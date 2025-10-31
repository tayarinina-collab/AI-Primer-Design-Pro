# -*- coding: utf-8 -*-
"""
AI Primer Design Pro – Advanced (Geneious-Pro Style)
Version 3.5 – mit KI-Parameteranalyse & robustem Primer3-Fallback
"""
import os
import io
import textwrap
import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio.Seq import Seq
from Bio import SeqIO

# Primer3-Integration
try:
    import primer3
    P3_OK = True
except Exception:
    P3_OK = False


# ========================== Hilfsfunktionen ===============================

def revcomp(seq: str) -> str:
    return str(Seq(seq).reverse_complement())


def gc_percent(seq: str) -> float:
    seq = seq.upper().replace("U", "T")
    return round(100 * (seq.count("G") + seq.count("C")) / len(seq), 2) if seq else 0.0


def longest_homopolymer(seq: str) -> int:
    max_run = run = 1
    for i in range(1, len(seq)):
        run = run + 1 if seq[i] == seq[i-1] else 1
        max_run = max(max_run, run)
    return max_run


def seq_stats_dna(seq: str) -> dict:
    seq = seq.upper().replace("U", "T")
    gc = gc_percent(seq)
    return {"len": len(seq), "gc": gc, "nA": seq.count("A"), "nC": seq.count("C"),
            "nG": seq.count("G"), "nT": seq.count("T")}


def dG_hairpin(seq: str) -> float:
    r = primer3.calcHairpin(seq)
    return r.dg if r.structure_found else 0.0


def dG_homodimer(seq: str) -> float:
    r = primer3.calcHomodimer(seq)
    return r.dg if r.structure_found else 0.0


def dG_heterodimer(a: str, b: str) -> float:
    r = primer3.calcHeterodimer(a, b)
    return r.dg if r.structure_found else 0.0


def simple_score(tm_l, tm_r, gc_l, gc_r, dg_self, dg_cross, homopoly, tm_target):
    tm_gap = abs(tm_l - tm_r)
    tm_mid = (tm_l + tm_r) / 2
    tm_pen = abs(tm_mid - tm_target)
    gc_pen = (abs(gc_l - 50) + abs(gc_r - 50)) / 2
    dim_pen = max(0, -min(dg_self, dg_cross))
    homo_pen = max(0, homopoly - 5) * 2
    score = 1 / (1 + 0.05 * tm_pen + 0.05 * tm_gap + 0.03 * gc_pen + 0.02 * dim_pen + 0.05 * homo_pen)
    return round(score, 3)


# ========================== KI-Parameteranalyse ===============================

def ki_param_vorschlaege(seq: str, args: dict, gc_clamp: bool, max_homopoly: int):
    """Offline-KI: gibt Vorschläge, wenn primer3 keine Primer findet"""
    stx = seq_stats_dna(seq)
    reasons = []
    sugg = args.copy()

    # Produktgröße (niemals < 60 bp)
    lo, hi = sugg["PRIMER_PRODUCT_SIZE_RANGE"][0]
    if stx["len"] < 150:
        lo2 = 60
        hi2 = max(100, min(400, stx["len"] - 20))
        sugg["PRIMER_PRODUCT_SIZE_RANGE"] = [[lo2, hi2]]
        reasons.append(f"Sequenz ist kurz ({stx['len']} bp) → Produktgröße angepasst auf {lo2}–{hi2} bp.")
    else:
        lo2 = min(80, lo)
        hi2 = max(600, hi)
        sugg["PRIMER_PRODUCT_SIZE_RANGE"] = [[lo2, hi2]]
        reasons.append(f"Produktgröße verbreitert auf {lo2}–{hi2} bp.")

    # GC%
    seq_gc = stx["gc"]
    gmin, gmax = sugg.get("PRIMER_MIN_GC", 40), sugg.get("PRIMER_MAX_GC", 60)
    if seq_gc < gmin or seq_gc > gmax or (gmax - gmin) < 20:
        span = 15
        new_min = int(max(30, seq_gc - span))
        new_max = int(min(75, seq_gc + span))
        sugg["PRIMER_MIN_GC"], sugg["PRIMER_MAX_GC"] = new_min, new_max
        reasons.append(f"GC-Bereich an Sequenz-GC ({seq_gc:.1f}%) angepasst → {new_min}–{new_max}%.")

    # Tm-Bereich
    tmin, tmax = sugg["PRIMER_MIN_TM"], sugg["PRIMER_MAX_TM"]
    if (tmax - tmin) < 10:
        tmid = (tmin + tmax) / 2
        sugg["PRIMER_MIN_TM"] = max(48, int(tmid - 7))
        sugg["PRIMER_MAX_TM"] = min(72, int(tmid + 7))
        reasons.append(f"Tm-Bereich verbreitert auf {sugg['PRIMER_MIN_TM']}–{sugg['PRIMER_MAX_TM']}°C.")

    # Primerlänge
    pmin, popt, pmax = sugg["PRIMER_MIN_SIZE"], sugg["PRIMER_OPT_SIZE"], sugg["PRIMER_MAX_SIZE"]
    sugg["PRIMER_MIN_SIZE"], sugg["PRIMER_OPT_SIZE"], sugg["PRIMER_MAX_SIZE"] = min(16, pmin), popt, max(32, pmax)
    reasons.append(f"Primerlängen erlaubt: {sugg['PRIMER_MIN_SIZE']},{sugg['PRIMER_OPT_SIZE']},{sugg['PRIMER_MAX_SIZE']} bp.")

    # GC-Clamp
    apply_gc_clamp = False if gc_clamp else gc_clamp
    if gc_clamp:
        reasons.append("3′-GC-Clamp temporär deaktiviert, um mehr Kandidaten zuzulassen.")

    # Kandidatenzahl erhöhen
    sugg["PRIMER_NUM_RETURN"] = max(36, int(sugg.get("PRIMER_NUM_RETURN", 24)) + 12)
    reasons.append(f"Kandidatenzahl auf {sugg['PRIMER_NUM_RETURN']} erhöht.")

    return sugg, apply_gc_clamp, reasons, stx


# ========================== Hauptfunktion ===============================

def run_primer_design_advanced():
    st.title("🧪 Primer Design – Advanced (Geneious-Pro Style)")
    st.caption("Automatisches Primer-Design, KI-Parameterhilfe, qPCR-Probe, Batch-Modus")

    # Eingabe
    seq_input = st.text_area("DNA-Sequenz (5'→3') eingeben oder FASTA laden:", height=150)
    if not seq_input.strip():
        st.info("Bitte Sequenz eingeben, um Primer zu entwerfen.")
        return

    seq = seq_input.upper().replace("U", "T").replace("\n", "")

    # Parameter
    st.subheader("⚙️ Design-Parameter")
    col1, col2, col3 = st.columns(3)
    with col1:
        primer_len = st.text_input("Primerlängen (min,opt,max)", "18,20,25")
        pmin, popt, pmax = [int(x) for x in primer_len.split(",")]
        prod_min, prod_max = st.slider("Produktgröße (bp)", 60, 1200, (80, 400))
    with col2:
        tm_min, tm_max = st.slider("Tm-Bereich (°C)", 48, 75, (58, 62))
        gc_min, gc_max = st.slider("GC-Bereich (%)", 30, 80, (40, 60))
    with col3:
        gc_clamp = st.checkbox("3′-GC-Clamp bevorzugen", True)
        max_homopoly = st.slider("Max. Homopolymer-Länge", 3, 8, 5)
        n_return = st.number_input("Anzahl zurückgegebener Primerpaare", 1, 50, 24)

    # Aktion
    if st.button("🚀 Primer entwerfen"):
        if not P3_OK:
            st.error("❌ primer3-py ist nicht installiert. Bitte hinzufügen.")
            return

        args = {
            "PRIMER_OPT_SIZE": popt,
            "PRIMER_MIN_SIZE": pmin,
            "PRIMER_MAX_SIZE": pmax,
            "PRIMER_MIN_TM": tm_min,
            "PRIMER_MAX_TM": tm_max,
            "PRIMER_MIN_GC": gc_min,
            "PRIMER_MAX_GC": gc_max,
            "PRIMER_PRODUCT_SIZE_RANGE": [list((prod_min, prod_max))],
            "PRIMER_NUM_RETURN": n_return,
            "PRIMER_EXPLAIN_FLAG": 1
        }

        design = primer3.bindings.designPrimers(
            {"SEQUENCE_ID": "target", "SEQUENCE_TEMPLATE": seq},
            args
        )

        total = design.get("PRIMER_PAIR_NUM_RETURNED", 0)

        if total == 0:
            st.warning("Keine Primer gefunden. KI-Assistent schlägt folgende Anpassungen vor:")
            sugg, gc_temp, reasons, stats = ki_param_vorschlaege(seq, args, gc_clamp, max_homopoly)

            df = pd.DataFrame([
                {"Parameter": "Tm (°C)", "Aktuell": f"{tm_min}–{tm_max}", "Vorschlag": f"{sugg['PRIMER_MIN_TM']}–{sugg['PRIMER_MAX_TM']}"},
                {"Parameter": "GC (%)", "Aktuell": f"{gc_min}–{gc_max}", "Vorschlag": f"{sugg['PRIMER_MIN_GC']}–{sugg['PRIMER_MAX_GC']}"},
                {"Parameter": "Primerlänge", "Aktuell": f"{pmin},{popt},{pmax}", "Vorschlag": f"{sugg['PRIMER_MIN_SIZE']},{sugg['PRIMER_OPT_SIZE']},{sugg['PRIMER_MAX_SIZE']}"},
                {"Parameter": "Produktgröße", "Aktuell": f"{prod_min}–{prod_max}", "Vorschlag": f"{sugg['PRIMER_PRODUCT_SIZE_RANGE'][0][0]}–{sugg['PRIMER_PRODUCT_SIZE_RANGE'][0][1]}"},
                {"Parameter": "Kandidatanzahl", "Aktuell": n_return, "Vorschlag": sugg['PRIMER_NUM_RETURN']}
            ])
            st.dataframe(df, use_container_width=True)
            st.markdown("**Begründungen:**")
            for r in reasons:
                st.write("• " + r)
            return

        # Ergebnisse anzeigen
        result = []
        for i in range(total):
            lp = design[f"PRIMER_LEFT_{i}_SEQUENCE"]
            rp = design[f"PRIMER_RIGHT_{i}_SEQUENCE"]
            tm_l, tm_r = design[f"PRIMER_LEFT_{i}_TM"], design[f"PRIMER_RIGHT_{i}_TM"]
            gc_l, gc_r = gc_percent(lp), gc_percent(rp)
            dg_self, dg_cross = dG_homodimer(lp), dG_heterodimer(lp, rp)
            homo = max(longest_homopolymer(lp), longest_homopolymer(rp))
            prod = design[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"]
            score = simple_score(tm_l, tm_r, gc_l, gc_r, dg_self, dg_cross, homo, (tm_min + tm_max) / 2)

            result.append({
                "Index": i + 1,
                "Left": lp,
                "Right": rp,
                "Tm (°C)": round(np.mean([tm_l, tm_r]), 1),
                "GC%": round(np.mean([gc_l, gc_r]), 1),
                "ΔG (kcal/mol)": round(min(dg_self, dg_cross), 2),
                "Homopolymer": homo,
                "Amplicon Size (bp)": prod,
                "Score": score
            })

        df = pd.DataFrame(result)
        st.success(f"{len(df)} Primerpaare erfolgreich generiert ✅")
        st.dataframe(df, use_container_width=True)

        # Visualisierung
        st.subheader("🧬 Visualisierung")
        fig, ax = plt.subplots(figsize=(6, 2))
        ax.scatter(df["GC%"], df["Tm (°C)"], c=df["Score"], cmap="viridis", s=80)
        ax.set_xlabel("GC (%)")
        ax.set_ylabel("Tm (°C)")
        ax.set_title("Primer-Thermodynamik-Profil")
        st.pyplot(fig)

        # Download
        st.download_button("⬇️ Exportiere Ergebnisse als CSV", df.to_csv(index=False).encode("utf-8"),
                           file_name="primer_results.csv", mime="text/csv")


# ========================== Ende Modul ===============================
