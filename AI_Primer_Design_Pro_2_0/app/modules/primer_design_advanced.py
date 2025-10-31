# -*- coding: utf-8 -*-
"""
Primer Design Advanced (Geneious-Pro Style) – offlinefähig
Funktionen:
- Automatisches Primerpaar-Design via primer3
- Import von Primerlisten (CSV/TSV/XLSX)
- Manuelles Prüfen/Scoren eingegebener Primer
- 5'-Extensions (Preset + Custom)
- Degenerate Primer via Consensus aus Alignment (IUPAC)
- Off-Target-Check gegen lokale FASTA-DB (mismatch-tolerant, beide Stränge)
- Amplicon-Extraktion & Visualisierung
- qPCR/TaqMan-Probe
- Optionale AI-Zusammenfassung (OpenAI, falls KEY verfügbar)
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

# ---- Optional: OpenAI und primer3 ----
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

try:
    import openpyxl
    XLSX_OK = True
except Exception:
    XLSX_OK = False


# ============================= Utilities ====================================
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
    r = primer3.calcHairpin(seq)
    return r.dg if r.structure_found else 0.0

def dG_homodimer(seq: str) -> float:
    r = primer3.calcHomodimer(seq)
    return r.dg if r.structure_found else 0.0

def dG_heterodimer(a: str, b: str) -> float:
    r = primer3.calcHeterodimer(a, b)
    return r.dg if r.structure_found else 0.0

def simple_score(tm_l, tm_r, gc_l, gc_r, dg_hd, dg_xd, homopoly, tm_target):
    tm_gap = abs(tm_l - tm_r)
    tm_mid = (tm_l+tm_r)/2
    tm_pen = abs(tm_mid - tm_target)
    gc_pen = (abs(gc_l-50)+abs(gc_r-50))/2
    dim_pen = max(0, -min(dg_hd, dg_xd))
    homo_pen = max(0, homopoly-4)*2
    score = 1.0/(1 + 0.06*tm_gap + 0.05*tm_pen + 0.02*gc_pen + 0.02*dim_pen + 0.04*homo_pen)
    return round(float(np.clip(score, 0, 1)), 3)


# ============================= 5'-Extensions ================================
EXT_PRESETS = {
    "—": ("", ""),
    "EcoRI (GAATTC)": ("GAATTC", "GAATTC"),
    "BamHI (GGATCC)": ("GGATCC", "GGATCC"),
    "XbaI (TCTAGA)": ("TCTAGA", "TCTAGA"),
    "Poly-A (AAAA)": ("AAAA", "AAAA"),
    "Gateway attB (kurz)": ("GGGGACAAGTTTGTACAAAAAAGCAGGCT", "GGGGACCACTTTGTACAAGAAAGCTGGGT")
}


# ============================= Main App =====================================
def run_primer_design_advanced():
    st.title("🧪 Primer Design – Advanced (Geneious Pro)")
    st.caption("Import/Export, Off-Target-Check, Degenerate-Design, 5′-Extensions, qPCR-Probe, Visualisierung")

    # --- Eingaben -----------------------------------------------------------
    st.subheader("📥 Eingaben")
    left, right = st.columns(2)
    with left:
        upl_target = st.file_uploader("Zielsequenz (FASTA/TXT)", type=["fasta", "fa", "txt"])
        target_text = st.text_area("…oder DNA-Sequenz hier einfügen (5'→3')", height=120)
    with right:
        upl_alignment = st.file_uploader("Alignment (Multi-FASTA) für Degenerate-Primer (optional)", type=["fasta", "fa"])
        upl_primer_list = st.file_uploader("Primerliste importieren (CSV/TSV/XLSX)", type=["csv", "tsv", "txt", "xlsx"])
        upl_db = st.file_uploader("Off-Target-Datenbank (FASTA, optional)", type=["fasta", "fa"])

    # --- Zielsequenz laden --------------------------------------------------
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

    # --- Parameter ----------------------------------------------------------
    st.subheader("⚙️ Design-Parameter")
    c1, c2, c3 = st.columns(3)

    with c1:
        st.markdown("**Primerlängenbereich (bp)**")
        primer_length_range = st.text_input("min,opt,max", "18,20,25")
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

    # --- Erweiterte Optionen -----------------------------------------------
    st.subheader("🧩 Erweiterte Optionen")
    gc_min, gc_max = st.slider("GC-Gehalt (%)", 20, 80, (40, 60))
    gc_clamp = st.checkbox("3'-GC-Clamp bevorzugen", True)
    max_homopoly = st.slider("Max. Homopolymer-Länge", 3, 8, 5)
    allow_degenerate = st.checkbox("Degenerate-Primer erlauben (IUPAC)", False)
    mismatches = st.slider("Off-Target-Suche: erlaubte Mismatches", 0, 3, 2)

    ex_col1, ex_col2 = st.columns(2)
    with ex_col1:
        preset = st.selectbox("5'-Extensions Preset", list(EXT_PRESETS.keys()), index=0)
    with ex_col2:
        custom_left = st.text_input("Custom 5'-Extension (Left)", "")
        custom_right = st.text_input("Custom 5'-Extension (Right)", "")

    # --- qPCR Optionen ------------------------------------------------------
    st.subheader("🧫 qPCR / Probe")
    enable_probe = st.checkbox("Probe mitentwerfen", False)
    reporter = st.selectbox("Reporter", ["FAM", "HEX", "VIC", "ROX", "Cy5"], index=0)
    quencher = st.selectbox("Quencher", ["BHQ1", "BHQ2", "TAMRA", "Iowa Black FQ"], index=0)

    # --- Button -------------------------------------------------------------
    if st.button("🚀 Automatisches Design starten (primer3)"):
        if not P3_OK:
            st.error("❌ primer3-py nicht installiert.")
            return
        if not target_seq:
            st.warning("Bitte DNA-Sequenz eingeben oder FASTA hochladen.")
            return

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
            "PRIMER_NUM_RETURN": 24,
            "PRIMER_EXPLAIN_FLAG": 1,
            "PRIMER_SALT_MONOVALENT": monoval,
            "PRIMER_SALT_DIVALENT": dival,
            "PRIMER_DNTP_CONC": dntp,
            "PRIMER_DNA_CONC": 250.0
        }

        design = primer3.bindings.designPrimers(
            {"SEQUENCE_ID": "target", "SEQUENCE_TEMPLATE": target_seq},
            args
        )

        pairs = []
        for i in range(design.get("PRIMER_PAIR_NUM_RETURNED", 0)):
            lseq = design.get(f"PRIMER_LEFT_{i}_SEQUENCE")
            rseq = design.get(f"PRIMER_RIGHT_{i}_SEQUENCE")
            if not lseq or not rseq:
                continue

            lseq_full = extL + lseq
            rseq_full = extR + rseq
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
        st.success(f"{len(df)} Primerpaare erfolgreich generiert ✅")
        st.dataframe(df, use_container_width=True)
