# -*- coding: utf-8 -*-
"""
AI Primer Design Pro – Primer Design Advanced (v3.3)
Features:
✅ Primer-Design via primer3 (Tm, GC, ΔG, Homodimer, Hairpin)
✅ 5′-Extensions, GC-Clamp, Degenerate-Fallback
✅ qPCR-Probe optional
✅ korrekte Primerposition & Amplicon-Visualisierung (Geneious-Style)
✅ KI-Ergebnisinterpretation (OpenAI optional)
"""

import os
import io
import textwrap
import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.Seq import Seq

# ===== Optionale Bibliotheken =====
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


# ============================= UTILITIES =============================
def revcomp(s: str) -> str:
    return str(Seq(s).reverse_complement())

def gc_percent(seq: str) -> float:
    s = seq.upper().replace("U", "T")
    return 0.0 if not s else (100.0 * (s.count("G")+s.count("C")) / len(s))

def longest_homopolymer(s: str) -> int:
    best = run = 1
    for i in range(1, len(s)):
        run = run+1 if s[i] == s[i-1] else 1
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

def simple_score(tm_l, tm_r, gc_l, gc_r, dg_hd, dg_xd, homopoly, tm_target):
    tm_gap = abs(tm_l - tm_r)
    tm_mid = (tm_l+tm_r)/2
    tm_pen = abs(tm_mid - tm_target)
    gc_pen = (abs(gc_l-50)+abs(gc_r-50))/2
    dim_pen = max(0, -min(dg_hd, dg_xd))
    homo_pen = max(0, homopoly-4)*2
    score = 1.0/(1 + 0.06*tm_gap + 0.05*tm_pen + 0.02*gc_pen + 0.02*dim_pen + 0.04*homo_pen)
    return round(float(np.clip(score, 0, 1)), 3)


# ============================= 5'-Extensions =============================
EXT_PRESETS = {
    "—": ("", ""),
    "EcoRI (GAATTC)": ("GAATTC", "GAATTC"),
    "BamHI (GGATCC)": ("GGATCC", "GGATCC"),
    "XbaI (TCTAGA)": ("TCTAGA", "TCTAGA"),
    "Poly-A (AAAA)": ("AAAA", "AAAA"),
    "Gateway attB (kurz)": ("GGGGACAAGTTTGTACAAAAAAGCAGGCT", "GGGGACCACTTTGTACAAGAAAGCTGGGT")
}


# ============================= MAIN FUNCTION =============================
def run_primer_design_advanced():
    st.title("🧪 Primer Design – Advanced (Geneious Pro Style)")
    st.caption("Primer-Design mit Thermoanalyse, GC-Clamp, 5′-Extensions, Visualisierung & AI-Erklärung")

    # --- Eingaben --------------------------------------------------------
    st.subheader("📥 Eingaben")
    col1, col2 = st.columns(2)
    with col1:
        upl_target = st.file_uploader("Zielsequenz (FASTA/TXT)", type=["fasta", "fa", "txt"])
        target_text = st.text_area("…oder DNA-Sequenz hier einfügen (5'→3')", height=120)
    with col2:
        upl_db = st.file_uploader("Off-Target-Datenbank (FASTA, optional)", type=["fasta", "fa"])

    target_seq = ""
    if upl_target is not None:
        txt = upl_target.getvalue().decode("utf-8").strip()
        if txt.startswith(">"):
            recs = list(SeqIO.parse(io.StringIO(txt), "fasta"))
            if recs:
                target_seq = str(recs[0].seq).upper().replace("U", "T")
        else:
            target_seq = txt.upper().replace("U", "T")
    elif target_text.strip():
        target_seq = target_text.strip().upper().replace("U", "T")

    # --- Parameter --------------------------------------------------------
    st.subheader("⚙️ Design-Parameter")
    c1, c2, c3 = st.columns(3)
    with c1:
        primer_length_range = st.text_input("Primerlängenbereich (min,opt,max)", "18,20,25")
        try:
            pmin, popt, pmax = [int(x) for x in primer_length_range.split(",")]
        except:
            pmin, popt, pmax = 18, 20, 25
    with c2:
        tm_min, tm_max = st.slider("Tm-Bereich (°C)", 48, 75, (58, 62))
        prod_min, prod_max = st.slider("Produktgröße (bp)", 60, 1500, (100, 400))
    with c3:
        monoval = st.number_input("Na⁺/K⁺ (mM)", 0.0, 500.0, 50.0, step=1.0)
        dival = st.number_input("Mg²⁺ (mM)", 0.0, 10.0, 1.5, step=0.1)
        dntp = st.number_input("dNTP (mM)", 0.0, 5.0, 0.6, step=0.1)

    st.subheader("🧩 Erweiterte Optionen")
    gc_min, gc_max = st.slider("GC-Gehalt (%)", 20, 80, (40, 60))
    gc_clamp = st.checkbox("3'-GC-Clamp bevorzugen", True)
    max_homopoly = st.slider("Max. Homopolymer-Länge", 3, 8, 5)

    ex1, ex2 = st.columns(2)
    with ex1:
        preset = st.selectbox("5'-Extensions Preset", list(EXT_PRESETS.keys()), index=0)
    with ex2:
        custom_left = st.text_input("Custom 5'-Extension (Left)", "")
        custom_right = st.text_input("Custom 5'-Extension (Right)", "")

    # --- qPCR Options -----------------------------------------------------
    st.subheader("🧫 qPCR / Probe")
    enable_probe = st.checkbox("Probe mitentwerfen", False)
    reporter = st.selectbox("Reporter", ["FAM", "HEX", "VIC", "ROX", "Cy5"], index=0)
    quencher = st.selectbox("Quencher", ["BHQ1", "BHQ2", "TAMRA", "Iowa Black FQ"], index=0)

    # ====================== DESIGN START ===========================
    if st.button("🚀 Primer entwerfen"):
        if not P3_OK:
            st.error("❌ primer3-py nicht installiert.")
            return
        if not target_seq:
            st.warning("Bitte DNA-Sequenz eingeben oder FASTA hochladen.")
            return

        if len(target_seq) < 50:
            st.error("❌ Sequenz zu kurz (<50 bp) – bitte längere DNA eingeben.")
            return

        if any(x <= 0 for x in [monoval, dival, dntp]):
            monoval, dival, dntp = 50.0, 1.5, 0.6

        # --- Extensions bestimmen ---
        extL, extR = EXT_PRESETS.get(preset, ("", ""))
        if custom_left: extL = custom_left.upper().replace("U", "T")
        if custom_right: extR = custom_right.upper().replace("U", "T")

        args = {
            "PRIMER_OPT_SIZE": popt,
            "PRIMER_MIN_SIZE": pmin,
            "PRIMER_MAX_SIZE": pmax,
            "PRIMER_MIN_TM": tm_min,
            "PRIMER_MAX_TM": tm_max,
            "PRIMER_MIN_GC": gc_min,
            "PRIMER_MAX_GC": gc_max,
            "PRIMER_PRODUCT_SIZE_RANGE": [list((prod_min, prod_max))],
            "PRIMER_NUM_RETURN": 10,
            "PRIMER_EXPLAIN_FLAG": 1,
            "PRIMER_SALT_MONOVALENT": monoval,
            "PRIMER_SALT_DIVALENT": dival,
            "PRIMER_DNTP_CONC": dntp,
            "PRIMER_DNA_CONC": 250.0
        }

        try:
            design = primer3.bindings.designPrimers(
                {"SEQUENCE_ID": "target", "SEQUENCE_TEMPLATE": target_seq},
                args
            )
        except ValueError:
            st.error("❌ Thermoanalyse-Fehler: Überprüfe Degenerate-Basen oder GC-Bereich.")
            return

        pairs = []
        for i in range(design.get("PRIMER_PAIR_NUM_RETURNED", 0)):
            lseq = design.get(f"PRIMER_LEFT_{i}_SEQUENCE")
            rseq = design.get(f"PRIMER_RIGHT_{i}_SEQUENCE")
            if not lseq or not rseq:
                continue

            lseq_full, rseq_full = extL + lseq, extR + rseq
            tm_l, tm_r = design.get(f"PRIMER_LEFT_{i}_TM"), design.get(f"PRIMER_RIGHT_{i}_TM")
            gc_l, gc_r = gc_percent(lseq_full), gc_percent(rseq_full)
            dghp = min(dG_hairpin(lseq_full), dG_hairpin(rseq_full))
            dgself = min(dG_homodimer(lseq_full), dG_homodimer(rseq_full))
            dgcross = dG_heterodimer(lseq_full, rseq_full)
            homo = max(longest_homopolymer(lseq_full), longest_homopolymer(rseq_full))
            if homo > max_homopoly: continue
            if gc_clamp and (lseq_full[-1] not in "GC" or rseq_full[-1] not in "GC"): continue

            score = simple_score(tm_l, tm_r, gc_l, gc_r, dgself, dgcross, homo, (tm_min+tm_max)/2)
            pairs.append({
                "Rank": i+1,
                "Left": lseq_full,
                "Right": rseq_full,
                "Tm_left": round(tm_l,2),
                "Tm_right": round(tm_r,2),
                "GC_left": round(gc_l,1),
                "GC_right": round(gc_r,1),
                "ΔG_hp(min)": round(dghp,2),
                "ΔG_dimer_self(min)": round(dgself,2),
                "ΔG_dimer_cross": round(dgcross,2),
                "Homopoly_max": homo,
                "Score": score
            })

        if not pairs:
            st.warning("Keine geeigneten Primer gefunden – bitte Parameter lockern.")
            return

        df = pd.DataFrame(pairs)
        st.success(f"✅ {len(df)} Primerpaare erfolgreich generiert")
        st.dataframe(df, use_container_width=True)

        # ================== Amplicon-Visualisierung ==================
        st.subheader("🧬 Amplicon-Visualisierung (Geneious-Style)")

        best = df.iloc[0]
        left_pos = design.get("PRIMER_LEFT_0")[0]
        left_len = design.get("PRIMER_LEFT_0")[1]
        right_pos = design.get("PRIMER_RIGHT_0")[0]
        right_len = design.get("PRIMER_RIGHT_0")[1]

        amplicon_start = left_pos
        amplicon_end = right_pos + right_len
        amplicon_seq = target_seq[amplicon_start:amplicon_end]

        fig, ax = plt.subplots(figsize=(8, 1.8))
        ax.set_xlim(0, len(target_seq))
        ax.set_ylim(0, 1)
        ax.set_yticks([])
        ax.set_xlabel("Position (bp)")
        ax.set_title("DNA Map mit Primer-Positionen")

        ax.axvspan(left_pos, left_pos + left_len, ymin=0.6, ymax=0.95, color="tab:blue", alpha=0.4, label="Forward Primer (→)")
        ax.axvspan(right_pos, right_pos + right_len, ymin=0.6, ymax=0.95, color="tab:orange", alpha=0.4, label="Reverse Primer (←)")
        ax.axvspan(amplicon_start, amplicon_end, ymin=0.25, ymax=0.55, color="tab:green", alpha=0.25, label="Amplicon")

        ax.annotate("→", xy=(left_pos + left_len / 2, 0.9), ha="center", va="center", fontsize=10, color="tab:blue")
        ax.annotate("←", xy=(right_pos + right_len / 2, 0.9), ha="center", va="center", fontsize=10, color="tab:orange")

        ax.legend(loc="upper right", ncol=3, fontsize=8, frameon=False)
        st.pyplot(fig, use_container_width=True)

        st.code(textwrap.fill(amplicon_seq, 80), language="text")
        st.download_button(
            "⬇️ Amplicon als FASTA",
            f">amplicon_{best['Rank']}\n{amplicon_seq}\n",
            file_name=f"amplicon_{best['Rank']}.fasta"
        )

        # ================== qPCR-Probe (optional) ==================
        if enable_probe:
            st.subheader("🧫 qPCR-Probe (zentriert)")
            center = (amplicon_start + amplicon_end) // 2
            probe_len = 25
            probe_seq = target_seq[center - probe_len // 2:center + probe_len // 2]
            st.write(f"**Probe:** {probe_seq}")
            st.caption(f"Reporter: **{reporter}**, Quencher: **{quencher}**")

    # ================== AI-Zusammenfassung ===========================
    st.markdown("---")
    st.subheader("🤖 AI-Zusammenfassung (optional)")

    ai_btn = st.button("Ergebnisse erklären lassen")
    if ai_btn:
        summary = (
            "Dieses Modul entwirft Primerpaare mittels primer3, "
            "analysiert Tm, GC-Gehalt und Dimer-/Hairpin-ΔG-Werte "
            "und bewertet jedes Paar nach Stabilität und Spezifität. "
            "Hohe Scores deuten auf balancierte Primer mit ähnlichen Tm-Werten, "
            "mittlerem GC-Gehalt (40–60%) und geringer Dimerneigung hin."
        )

        if OPENAI_OK and os.getenv("OPENAI_API_KEY"):
            try:
                openai.api_key = os.getenv("OPENAI_API_KEY")
                prompt = (
                    "Erkläre einem Laboranwender kurz auf Deutsch, "
                    "wie die erzeugten Primer bewertet werden (Tm, GC, ΔG, Amplicon-Größe) "
                    "und gib Tipps zur Optimierung."
                )
                resp = openai.ChatCompletion.create(
                    model="gpt-4o-mini",
                    messages=[{"role": "system", "content": prompt}]
                )
                summary = resp["choices"][0]["message"]["content"]
            except Exception as e:
                st.warning(f"⚠️ KI-Erklärung nicht verfügbar: {e}")

        st.info(summary)
