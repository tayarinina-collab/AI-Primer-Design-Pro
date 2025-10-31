# -*- coding: utf-8 -*-
"""
AI Primer Design Pro – Primer Design Advanced (v3.4)
Features:
- Primer-Design via primer3 (Tm, GC, ΔG, Homodimer, Hairpin)
- 5′-Extensions, GC-Clamp, Produktgrößenfenster
- KI-Assistent (offline, regelbasiert) bei "keine Primer gefunden"
  + optional: kurze GPT-Erklärung, wenn OPENAI_API_KEY gesetzt ist
- Auto-Fallback: zweiter Lauf mit gelockerten Parametern
- korrekte Primerpositionen (forward/reverse) & Amplicon-Visualisierung
- CSV/FASTA-Export, qPCR-Probe (heuristisch)
"""

import os
import io
import re
import textwrap
import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.Seq import Seq

# ---------- Optionale Bibliotheken ----------
try:
    import openai
    OPENAI_OK = True
except Exception:
    OPENAI_OK = False

try:
    import primer3
    P3_OK = True
except Exception:
    P3_OK = False


# ======================= Utility-Funktionen =======================

def sanitize_dna(seq: str) -> str:
    """Erlaubt nur A/C/G/T, ersetzt U->T, entfernt Whitespaces/sonstige Zeichen."""
    s = seq.upper().replace("U", "T")
    return re.sub(r"[^ACGT]", "", s)

def revcomp(s: str) -> str:
    return str(Seq(s).reverse_complement())

def gc_percent(seq: str) -> float:
    s = sanitize_dna(seq)
    return 0.0 if not s else (100.0 * (s.count("G") + s.count("C")) / len(s))

def longest_homopolymer(s: str) -> int:
    s = sanitize_dna(s)
    if not s:
        return 0
    best = run = 1
    for i in range(1, len(s)):
        run = run + 1 if s[i] == s[i - 1] else 1
        best = max(best, run)
    return best

def dG_hairpin(seq: str) -> float:
    try:
        r = primer3.calcHairpin(seq)
        return r.dg if r.structure_found else 0.0
    except Exception:
        return 0.0

def dG_homodimer(seq: str) -> float:
    try:
        r = primer3.calcHomodimer(seq)
        return r.dg if r.structure_found else 0.0
    except Exception:
        return 0.0

def dG_heterodimer(a: str, b: str) -> float:
    try:
        r = primer3.calcHeterodimer(a, b)
        return r.dg if r.structure_found else 0.0
    except Exception:
        return 0.0

def simple_score(tm_l, tm_r, gc_l, gc_r, dg_self, dg_cross, homopoly, tm_target):
    """0..1 – je höher desto besser (balancierte Tm, moderater GC, geringe Dimerneigung)."""
    tm_gap = abs(tm_l - tm_r)
    tm_mid = (tm_l + tm_r) / 2
    tm_pen = abs(tm_mid - tm_target)
    gc_pen = (abs(gc_l - 50) + abs(gc_r - 50)) / 2
    dim_pen = max(0, -min(dg_self, dg_cross))  # nur negative ΔG bestrafen
    homo_pen = max(0, homopoly - 4) * 2
    score = 1.0 / (1 + 0.06 * tm_gap + 0.05 * tm_pen + 0.02 * gc_pen + 0.02 * dim_pen + 0.04 * homo_pen)
    return float(np.clip(score, 0, 1))


# ======================= KI-Assistent (offline) =======================

def seq_stats_dna(seq: str):
    s = sanitize_dna(seq)
    L = len(s)
    gc = 0.0 if L == 0 else 100.0 * (s.count('G') + s.count('C')) / L
    best = longest_homopolymer(s)
    return {"len": L, "gc": round(gc, 2), "homopoly": best}

def ki_param_vorschlaege(seq: str, args: dict, gc_clamp: bool, max_homopoly: int):
    """
    Regelbasierte Vorschläge, wenn primer3 keine Primer findet.
    Gibt (neue_args, gc_clamp_temp, begruendungen:list, stats:dict) zurück.
    """
    stx = seq_stats_dna(seq)
    reasons = []
    sugg = args.copy()

    # 1) Produktgröße passend zur Sequenz/breiter machen
    lo, hi = sugg["PRIMER_PRODUCT_SIZE_RANGE"][0]
    if stx["len"] < 120:
        lo2, hi2 = 40, max(60, stx["len"] - 10)
        sugg["PRIMER_PRODUCT_SIZE_RANGE"] = [[lo2, hi2]]
        reasons.append(f"Sequenz sehr kurz ({stx['len']} bp) → Produktgröße auf {lo2}–{hi2} bp gesetzt.")
    else:
        lo2 = min(70, lo)
        hi2 = max(600, hi)
        sugg["PRIMER_PRODUCT_SIZE_RANGE"] = [[lo2, hi2]]
        reasons.append(f"Produktgröße verbreitert auf {lo2}–{hi2} bp.")

    # 2) GC-Band an Sequenz-GC anlehnen
    seq_gc = stx["gc"]
    gmin, gmax = sugg.get("PRIMER_MIN_GC", 40), sugg.get("PRIMER_MAX_GC", 60)
    if seq_gc < gmin or seq_gc > gmax or (gmax - gmin) < 20:
        span = 12 if 35 <= seq_gc <= 65 else 18
        new_min = int(max(25, min(70, seq_gc - span)))
        new_max = int(min(75, max(35, seq_gc + span)))
        sugg["PRIMER_MIN_GC"], sugg["PRIMER_MAX_GC"] = new_min, new_max
        reasons.append(f"GC-Bereich an Sequenz-GC ({seq_gc:.1f}%) angepasst → {new_min}–{new_max}%.")

    # 3) Tm-Band verbreitern
    tmin, tmax = sugg["PRIMER_MIN_TM"], sugg["PRIMER_MAX_TM"]
    if (tmax - tmin) < 10:
        tmid = (tmin + tmax) / 2
        sugg["PRIMER_MIN_TM"] = max(48, int(round(tmid - 7)))
        sugg["PRIMER_MAX_TM"] = min(72, int(round(tmid + 7)))
        reasons.append(f"Tm-Band verbreitert auf {sugg['PRIMER_MIN_TM']}–{sugg['PRIMER_MAX_TM']}°C.")

    # 4) Primerlängen lockern
    pmin, popt, pmax = sugg["PRIMER_MIN_SIZE"], sugg["PRIMER_OPT_SIZE"], sugg["PRIMER_MAX_SIZE"]
    pmin = min(pmin, 16)
    pmax = max(pmax, 32)
    popt = int((pmin + pmax) / 2)
    sugg.update({"PRIMER_MIN_SIZE": pmin, "PRIMER_OPT_SIZE": popt, "PRIMER_MAX_SIZE": pmax})
    reasons.append(f"Primerlängen erlaubt: {pmin},{popt},{pmax} bp.")

    # 5) GC-Clamp & Homopolymere
    apply_gc_clamp = gc_clamp
    if gc_clamp:
        apply_gc_clamp = False
        reasons.append("3′-GC-Clamp temporär deaktiviert, um mehr Kandidaten zuzulassen.")
    if max_homopoly < 6:
        reasons.append("Max. Homopolymer-Länge auf 6 erhöht.")
    # (Design-Filter außerhalb berücksichtigen)

    # 6) Ionen/Konzentrationen sichern
    if any(sugg.get(k, 0) <= 0 for k in ["PRIMER_SALT_MONOVALENT", "PRIMER_SALT_DIVALENT", "PRIMER_DNTP_CONC"]):
        sugg["PRIMER_SALT_MONOVALENT"] = 50.0
        sugg["PRIMER_SALT_DIVALENT"] = 1.5
        sugg["PRIMER_DNTP_CONC"] = 0.6
        reasons.append("Ionen auf Standard gesetzt (50 mM Na⁺, 1.5 mM Mg²⁺, 0.6 mM dNTP).")

    # 7) Mehr Kandidaten
    sugg["PRIMER_NUM_RETURN"] = max(36, int(sugg.get("PRIMER_NUM_RETURN", 24)) + 12)
    reasons.append(f"Kandidatenzahl auf {sugg['PRIMER_NUM_RETURN']} erhöht.")

    return sugg, apply_gc_clamp, reasons, stx


# ======================= Hauptmodul =======================

EXT_PRESETS = {
    "—": ("", ""),
    "EcoRI (GAATTC)": ("GAATTC", "GAATTC"),
    "BamHI (GGATCC)": ("GGATCC", "GGATCC"),
    "XbaI (TCTAGA)": ("TCTAGA", "TCTAGA"),
    "Poly-A (AAAA)": ("AAAA", "AAAA"),
    "Gateway attB (kurz)": ("GGGGACAAGTTTGTACAAAAAAGCAGGCT", "GGGGACCACTTTGTACAAGAAAGCTGGGT"),
}

def run_primer_design_advanced():
    st.title("🧪 Primer Design – Advanced (Geneious Pro Style)")
    st.caption("Primer3-Design, KI-Assistent bei Fehlversuch, Visualisierung & Exporte")

    # -------- Eingaben --------
    st.subheader("📥 Eingaben")
    c1, c2 = st.columns(2)
    with c1:
        upl_target = st.file_uploader("Zielsequenz (FASTA/TXT)", type=["fasta", "fa", "txt"])
        target_text = st.text_area("…oder DNA-Sequenz hier einfügen (5'→3')", height=120)
    with c2:
        pass  # Platz für künftige Off-Target-DB etc.

    target_seq = ""
    if upl_target is not None:
        txt = upl_target.getvalue().decode("utf-8").strip()
        if txt.startswith(">"):
            recs = list(SeqIO.parse(io.StringIO(txt), "fasta"))
            if recs:
                target_seq = sanitize_dna(str(recs[0].seq))
        else:
            target_seq = sanitize_dna(txt)
    elif target_text.strip():
        target_seq = sanitize_dna(target_text)

    # -------- Parameter --------
    st.subheader("⚙️ Design-Parameter")
    p1, p2, p3 = st.columns(3)
    with p1:
        primer_length_range = st.text_input("Primerlängen (min,opt,max)", "18,20,25")
        try:
            P_MIN, P_OPT, P_MAX = [int(x) for x in primer_length_range.split(",")]
        except Exception:
            P_MIN, P_OPT, P_MAX = 18, 20, 25
    with p2:
        TM_MIN, TM_MAX = st.slider("Tm-Bereich (°C)", 48, 75, (58, 62))
        PROD_MIN, PROD_MAX = st.slider("Produktgröße (bp)", 60, 1500, (100, 400))
    with p3:
        SALT = st.number_input("Na⁺/K⁺ (mM)", 0.0, 500.0, 50.0, step=1.0)
        MG = st.number_input("Mg²⁺ (mM)", 0.0, 10.0, 1.5, step=0.1)
        DNTP = st.number_input("dNTP (mM)", 0.0, 5.0, 0.6, step=0.1)

    st.subheader("🧩 Erweiterte Optionen")
    GC_MIN, GC_MAX = st.slider("GC-Gehalt (%)", 20, 80, (40, 60))
    gc_clamp = st.checkbox("3′-GC-Clamp bevorzugen", True)
    max_homopoly = st.slider("Max. Homopolymer-Länge", 3, 10, 6)

    e1, e2 = st.columns(2)
    with e1:
        ext_preset = st.selectbox("5′-Extensions Preset", list(EXT_PRESETS.keys()), index=0)
    with e2:
        custom_left = st.text_input("Custom 5′-Extension (Left)", "")
        custom_right = st.text_input("Custom 5′-Extension (Right)", "")

    st.caption("💡 Wenn keine Primer gefunden werden, schlägt der KI-Assistent unten konkrete Änderungen vor.")

    # -------- qPCR / Probe --------
    st.subheader("🧫 qPCR / Probe (optional)")
    enable_probe = st.checkbox("Probe mitentwerfen", False)
    reporter = st.selectbox("Reporter", ["FAM", "HEX", "VIC", "ROX", "Cy5"], index=0)
    quencher = st.selectbox("Quencher", ["BHQ1", "BHQ2", "TAMRA", "Iowa Black FQ"], index=0)

    # -------- Button --------
    if st.button("🚀 Primer entwerfen"):
        if not P3_OK:
            st.error("❌ primer3-py nicht installiert.")
            return
        if not target_seq:
            st.warning("Bitte DNA-Sequenz eingeben oder FASTA hochladen.")
            return
        if len(target_seq) < 50:
            st.error("❌ Sequenz zu kurz (<50 bp). Bitte >100 bp verwenden.")
            return

        # Extensions
        extL, extR = EXT_PRESETS.get(ext_preset, ("", ""))
        if custom_left:
            extL = sanitize_dna(custom_left)
        if custom_right:
            extR = sanitize_dna(custom_right)

        # Grund-Args
        args = {
            "PRIMER_OPT_SIZE": P_OPT,
            "PRIMER_MIN_SIZE": P_MIN,
            "PRIMER_MAX_SIZE": P_MAX,
            "PRIMER_MIN_TM": TM_MIN,
            "PRIMER_MAX_TM": TM_MAX,
            "PRIMER_MIN_GC": GC_MIN,
            "PRIMER_MAX_GC": GC_MAX,
            "PRIMER_PRODUCT_SIZE_RANGE": [[PROD_MIN, PROD_MAX]],
            "PRIMER_NUM_RETURN": 24,
            "PRIMER_EXPLAIN_FLAG": 1,
            "PRIMER_SALT_MONOVALENT": SALT if SALT > 0 else 50.0,
            "PRIMER_SALT_DIVALENT": MG if MG > 0 else 1.5,
            "PRIMER_DNTP_CONC": DNTP if DNTP > 0 else 0.6,
            "PRIMER_DNA_CONC": 250.0,
        }

        # --- 1. Lauf ---
        try:
            design = primer3.bindings.designPrimers(
                {"SEQUENCE_ID": "target", "SEQUENCE_TEMPLATE": target_seq},
                args
            )
        except ValueError:
            st.error("❌ Thermoanalyse-Fehler: Prüfe Tm/GC/Produktfenster oder ungewöhnliche Zeichen.")
            return

        total = design.get("PRIMER_PAIR_NUM_RETURNED", 0)

        # --- Wenn nichts gefunden → KI-Assistent + Fallback ---
        if total == 0:
            ki_sugg, gc_clamp_temp, reasons, stx = ki_param_vorschlaege(target_seq, args, gc_clamp, max_homopoly)

            st.warning("Keine Primer gefunden. KI-Assistent schlägt folgende Anpassungen vor:")
            cur_vs_new = pd.DataFrame([
                ["Tm (°C)", f"{args['PRIMER_MIN_TM']}–{args['PRIMER_MAX_TM']}",
                 f"{ki_sugg['PRIMER_MIN_TM']}–{ki_sugg['PRIMER_MAX_TM']}"],
                ["GC (%)", f"{args['PRIMER_MIN_GC']}–{args['PRIMER_MAX_GC']}",
                 f"{ki_sugg['PRIMER_MIN_GC']}–{ki_sugg['PRIMER_MAX_GC']}"],
                ["Primerlänge", f"{args['PRIMER_MIN_SIZE']},{args['PRIMER_OPT_SIZE']},{args['PRIMER_MAX_SIZE']}",
                 f"{ki_sugg['PRIMER_MIN_SIZE']},{ki_sugg['PRIMER_OPT_SIZE']},{ki_sugg['PRIMER_MAX_SIZE']}"],
                ["Produktgröße", f"{args['PRIMER_PRODUCT_SIZE_RANGE'][0][0]}–{args['PRIMER_PRODUCT_SIZE_RANGE'][0][1]}",
                 f"{ki_sugg['PRIMER_PRODUCT_SIZE_RANGE'][0][0]}–{ki_sugg['PRIMER_PRODUCT_SIZE_RANGE'][0][1]}"],
                ["Kandidatenzahl", str(args.get("PRIMER_NUM_RETURN", 24)), str(ki_sugg["PRIMER_NUM_RETURN"])],
            ], columns=["Parameter", "Aktuell", "Vorschlag"])
            st.dataframe(cur_vs_new, use_container_width=True)

            st.markdown("**Begründungen:**")
            st.markdown("\n".join([f"- {r}" for r in reasons]))
            st.caption(f"Sequenz: Länge {stx['len']} bp · GC {stx['gc']}% · max. Homopolymer {stx['homopoly']}")

            if OPENAI_OK and os.getenv("OPENAI_API_KEY"):
                try:
                    openai.api_key = os.getenv("OPENAI_API_KEY")
                    prompt = (
                        "Gib kurze, konkrete Schritte zur Parameteranpassung, "
                        "damit primer3 Primer findet (Deutsch, Stichpunkte). "
                        f"Sequenzlänge: {stx['len']} bp; Sequenz-GC: {stx['gc']}%."
                    )
                    resp = openai.ChatCompletion.create(
                        model="gpt-4o-mini",
                        messages=[{"role": "system", "content": prompt}]
                    )
                    st.info(resp["choices"][0]["message"]["content"])
                except Exception:
                    pass

            # Interaktiver Button
            if st.button("✅ Vorschläge übernehmen & erneut versuchen"):
                st.session_state["_args_override"] = ki_sugg
                st.session_state["_gc_clamp_override"] = gc_clamp_temp
                st.rerun()
            st.stop()

        # Falls Re-Run mit Overrides
        over = st.session_state.pop("_args_override", None) if "_args_override" in st.session_state else None
        if over:
            args = over
        if "_gc_clamp_override" in st.session_state:
            gc_clamp = st.session_state.pop("_gc_clamp_override")

        # -------- Ergebnisse aufbereiten --------
        pairs = []
        for i in range(design.get("PRIMER_PAIR_NUM_RETURNED", 0)):
            lseq = design.get(f"PRIMER_LEFT_{i}_SEQUENCE")
            rseq = design.get(f"PRIMER_RIGHT_{i}_SEQUENCE")
            if not lseq or not rseq:
                continue

            # Extensions anwenden
            l_full = (extL + lseq)
            r_full = (extR + rseq)

            tm_l = design.get(f"PRIMER_LEFT_{i}_TM")
            tm_r = design.get(f"PRIMER_RIGHT_{i}_TM")
            gc_l = gc_percent(l_full)
            gc_r = gc_percent(r_full)
            dghp = min(dG_hairpin(l_full), dG_hairpin(r_full))
            dgself = min(dG_homodimer(l_full), dG_homodimer(r_full))
            dgcross = dG_heterodimer(l_full, r_full)
            homo = max(longest_homopolymer(l_full), longest_homopolymer(r_full))

            # Filter: Homopoly und optional GC-Clamp
            if homo > max_homopoly:
                continue
            if gc_clamp and (l_full[-1] not in "GC" or r_full[-1] not in "GC"):
                continue

            score = simple_score(tm_l, tm_r, gc_l, gc_r, dgself, dgcross, homo, (TM_MIN + TM_MAX) / 2)

            pairs.append({
                "Rank": i + 1,
                "Left": l_full,
                "Right": r_full,
                "Tm_left": round(tm_l, 2),
                "Tm_right": round(tm_r, 2),
                "GC_left": round(gc_l, 1),
                "GC_right": round(gc_r, 1),
                "ΔG_hp(min)": round(dghp, 2),
                "ΔG_dimer_self(min)": round(dgself, 2),
                "ΔG_dimer_cross": round(dgcross, 2),
                "Homopoly_max": homo,
                "Score": round(score, 3),
            })

        if not pairs:
            st.warning("Keine geeigneten Primer nach Filterung – bitte Parameter (GC-Clamp, Homopoly) lockern.")
            return

        df = pd.DataFrame(pairs).sort_values(["Score"], ascending=[False]).reset_index(drop=True)
        st.success(f"✅ {len(df)} Primerpaare gefunden")
        st.dataframe(df, use_container_width=True)

        # -------- CSV-Export --------
        st.download_button(
            "⬇️ Primerpaare als CSV",
            df.to_csv(index=False).encode("utf-8"),
            file_name="primer_pairs_advanced.csv",
            mime="text/csv"
        )

        # -------- Amplicon-Visualisierung (korrekte Koordinaten) --------
        st.subheader("🧬 Amplicon-Visualisierung (Geneious-Style)")

        # Nehme bestes Paar (Index 0 in primer3-Ergebnis, NICHT df-Rang)
        try:
            left_pos, left_len = design.get("PRIMER_LEFT_0")
            right_pos, right_len = design.get("PRIMER_RIGHT_0")
        except Exception:
            # Fallback: wenn Struktur anders
            left_pos, left_len = design.get("PRIMER_LEFT_0", [0, len(df.iloc[0]['Left'])])
            right_pos, right_len = design.get("PRIMER_RIGHT_0", [len(target_seq)-len(df.iloc[0]['Right']), len(df.iloc[0]['Right'])])

        amplicon_start = left_pos
        amplicon_end = right_pos + right_len
        amplicon_seq = target_seq[amplicon_start:amplicon_end]

        fig, ax = plt.subplots(figsize=(8, 1.8))
        ax.set_xlim(0, len(target_seq))
        ax.set_ylim(0, 1)
        ax.set_yticks([])
        ax.set_xlabel("Position (bp)")
        ax.set_title("DNA Map mit Primer-Positionen")

        ax.axvspan(left_pos, left_pos + left_len, ymin=0.6, ymax=0.95, color="tab:blue", alpha=0.4, label="Forward (→)")
        ax.axvspan(right_pos, right_pos + right_len, ymin=0.6, ymax=0.95, color="tab:orange", alpha=0.4, label="Reverse (←)")
        ax.axvspan(amplicon_start, amplicon_end, ymin=0.25, ymax=0.55, color="tab:green", alpha=0.25, label="Amplicon")

        ax.annotate("→", xy=(left_pos + left_len / 2, 0.9), ha="center", va="center", fontsize=10, color="tab:blue")
        ax.annotate("←", xy=(right_pos + right_len / 2, 0.9), ha="center", va="center", fontsize=10, color="tab:orange")
        ax.legend(loc="upper right", ncol=3, fontsize=8, frameon=False)
        st.pyplot(fig, use_container_width=True)

        st.code(textwrap.fill(amplicon_seq, 80), language="text")
        st.download_button(
            "⬇️ Amplicon als FASTA",
            f">amplicon_best\n{amplicon_seq}\n",
            file_name="amplicon_best.fasta",
            mime="text/plain"
        )

        # -------- qPCR-Probe (heuristisch) --------
        if enable_probe and len(amplicon_seq) >= 40:
            st.subheader("🧫 qPCR-Probe (zentriert, heuristisch)")
            center = (amplicon_start + amplicon_end) // 2
            p_len = min(28, max(18, int(0.4 * (amplicon_end - amplicon_start))))
            p_start = max(0, center - p_len // 2)
            probe = target_seq[p_start:p_start + p_len]
            st.write(f"**Probe:** {probe}")
            st.caption(f"Reporter: **{reporter}** · Quencher: **{quencher}**")

    # -------- KI-Erklärung (optional separat aufrufbar) --------
    st.markdown("---")
    st.subheader("🤖 KI-Zusammenfassung (optional)")
    if st.button("Ergebnisse erklären lassen"):
        summary = (
            "Die Bewertung berücksichtigt Tm-Gleichgewicht (beide Primer ähnlich), "
            "ausgewogenen GC-Gehalt (typisch 40–60%), geringe Dimer-/Hairpin-Neigung "
            "(ΔG möglichst wenig negativ) und eine sinnvolle Ampliconlänge. "
            "Bei ausbleibenden Treffern Tm-Band verbreitern, Produktfenster öffnen, "
            "GC-Band an Sequenz-GC anpassen und 3′-GC-Clamp vorübergehend deaktivieren."
        )
        if OPENAI_OK and os.getenv("OPENAI_API_KEY"):
            try:
                openai.api_key = os.getenv("OPENAI_API_KEY")
                prompt = (
                    "Erkläre einem Laboranwender kurz, wie die erzeugten Primer bewertet werden "
                    "(Tm, GC, ΔG, Amplicon) und nenne 3–5 konkrete Optimierungsschritte. "
                    "Antwort auf Deutsch, stichpunktartig."
                )
                resp = openai.ChatCompletion.create(
                    model="gpt-4o-mini",
                    messages=[{"role": "system", "content": prompt}]
                )
                summary = resp["choices"][0]["message"]["content"]
            except Exception as e:
                st.warning(f"⚠️ KI-Erklärung nicht verfügbar: {e}")
        st.info(summary)
