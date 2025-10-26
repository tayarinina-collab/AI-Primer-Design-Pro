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
import csv
import math
import json
import textwrap
import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ---- optionale AI (läuft nur, wenn OPENAI_API_KEY gesetzt ist) -------------
try:
    import openai
    OPENAI_OK = True
except Exception:
    OPENAI_OK = False

# ---- primer3 & biopython ---------------------------------------------------
try:
    import primer3
    P3_OK = True
except Exception:
    P3_OK = False

from Bio import SeqIO
from Bio.Seq import Seq

# ---- (empfohlen) Excel-Support, falls installiert --------------------------
try:
    import openpyxl  # noqa
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

# ---- Konsensus mit IUPAC aus Alignment (Multi-FASTA) -----------------------
IUPAC_MERGE = {
    frozenset("A"): "A", frozenset("C"): "C", frozenset("G"): "G", frozenset("T"): "T",
    frozenset("AG"): "R", frozenset("CT"): "Y", frozenset("GC"): "S", frozenset("AT"): "W",
    frozenset("GT"): "K", frozenset("AC"): "M",
    frozenset("CGT"): "B", frozenset("AGT"): "D", frozenset("ACT"): "H", frozenset("ACG"): "V",
    frozenset("ACGT"): "N"
}
def iupac_consensus_from_alignment(records):
    """records: list(SeqRecord) gleiche Länge (MSA)."""
    if not records:
        return ""
    L = len(records[0].seq)
    cons = []
    for i in range(L):
        bases = set()
        for r in records:
            b = r.seq[i].upper().replace("U", "T")
            if b in "-N":  # Lücken ignorieren
                continue
            if b not in "ACGT":
                b = "N"
            bases.add(b)
        if not bases:
            cons.append("N")
        else:
            cons.append(IUPAC_MERGE.get(frozenset(sorted(bases)), "N"))
    return "".join(cons)

# ---- Off-target Suche (mismatch-tolerant, beide Stränge) -------------------
def hamming(a: str, b: str) -> int:
    return sum(1 for x,y in zip(a,b) if x!=y)

def count_offtargets(primer: str, fasta_records, mismatches: int = 2) -> int:
    """Naiver Off-target-Scan über alle Sequenzen und beide Stränge."""
    primer = primer.upper()
    rc = revcomp(primer)
    k = len(primer)
    hits = 0
    for rec in fasta_records:
        s = str(rec.seq).upper()
        for strand in (s, revcomp(s)):
            if len(strand) < k:
                continue
            for i in range(0, len(strand)-k+1):
                if hamming(primer, strand[i:i+k]) <= mismatches:
                    hits += 1
                if hamming(rc, strand[i:i+k]) <= mismatches:
                    hits += 1
    return hits

# ---- 5'-Extensions presets --------------------------------------------------
EXT_PRESETS = {
    "—": ("",""),
    "EcoRI (GAATTC)": ("GAATTC", "GAATTC"),
    "BamHI (GGATCC)": ("GGATCC", "GGATCC"),
    "XbaI (TCTAGA)": ("TCTAGA", "TCTAGA"),
    "Poly-A (AAAA)": ("AAAA", "AAAA"),
    "Gateway attB (kurz)": ("GGGGACAAGTTTGTACAAAAAAGCAGGCT", "GGGGACCACTTTGTACAAGAAAGCTGGGT")
}

# ============================== UI ==========================================
def run_primer_design_advanced():
    st.title("🧪 Primer Design – Advanced (Geneious Pro)")
    st.caption("Import/Export, Off-Target-Check, Degenerate-Design, 5′-Extensions, qPCR-Probe, Visualisierung")

    # --- Eingaben: Zielsequenz / Alignment / Primerlisten / Off-target DB ----
    st.subheader("📥 Eingaben")
    left, right = st.columns(2)
    with left:
        upl_target = st.file_uploader("Zielsequenz (FASTA/TXT)", type=["fasta","fa","txt"])
        target_text = st.text_area("…oder DNA-Sequenz hier einfügen (5'→3')", height=120)
    with right:
        upl_alignment = st.file_uploader("Alignment (Multi-FASTA) für Degenerate-Primer (optional)", type=["fasta","fa"])
        upl_primer_list = st.file_uploader("Primerliste importieren (CSV/TSV/XLSX)", type=["csv","tsv","txt","xlsx"])
        upl_db = st.file_uploader("Off-Target-Datenbank (FASTA, optional)", type=["fasta","fa"])

    # Zielsequenz laden
    target_seq = ""
    if upl_target is not None:
        txt = upl_target.getvalue().decode("utf-8").strip()
        if txt.startswith(">"):
            recs = list(SeqIO.parse(io.StringIO(txt), "fasta"))
            if recs:
                target_seq = str(recs[0].seq).upper().replace("U","T")
        else:
            target_seq = txt.upper().replace("U","T")
    elif target_text.strip():
        target_seq = target_text.strip().upper().replace("U","T")

    # Alignment-Konsensus für Degenerate
    consensus = ""
    if upl_alignment is not None:
        try:
            recs = list(SeqIO.parse(io.StringIO(upl_alignment.getvalue().decode("utf-8")), "fasta"))
            consensus = iupac_consensus_from_alignment(recs)
        except Exception:
            st.warning("⚠️ Konnte Alignment nicht lesen (Multi-FASTA erwartet).")

    # Primerliste importieren
    imported_oligos = pd.DataFrame()
    if upl_primer_list is not None:
        try:
            if upl_primer_list.name.lower().endswith(".xlsx"):
                if XLSX_OK:
                    imported_oligos = pd.read_excel(upl_primer_list)
                else:
                    st.warning("⚠️ Für XLSX bitte `openpyxl` in requirements hinzufügen.")
            elif upl_primer_list.name.lower().endswith(".tsv"):
                imported_oligos = pd.read_csv(upl_primer_list, sep="\t")
            else:
                imported_oligos = pd.read_csv(upl_primer_list)
            st.success(f"🔬 {len(imported_oligos)} Primer aus Datei geladen.")
        except Exception as e:
            st.error(f"Fehler beim Import: {e}")

    # Off-target DB laden
    offtarget_db = []
    if upl_db is not None:
        try:
            offtarget_db = list(SeqIO.parse(io.StringIO(upl_db.getvalue().decode("utf-8")), "fasta"))
            st.info(f"🔎 Off-Target-DB: {len(offtarget_db)} Sequenzen geladen.")
        except Exception:
            st.warning("⚠️ Konnte Off-Target-FASTA nicht lesen.")

    # -------------------- Parameter ------------------------------------------
  st.subheader("⚙️ Design-Parameter")

c1, c2, c3 = st.columns(3)

with c1:
    st.markdown("**Primerlängenbereich (bp)**")
    primer_length_range = st.text_input("min,opt,max", "18,20,25")
    try:
        pmin, popt, pmax = [int(x) for x in primer_length_range.split(",")]
    except:
        st.warning("⚠️ Bitte drei Werte eingeben, z. B. 18,20,25")
        pmin, popt, pmax = 18, 20, 25  # Standardwerte als Fallback

with c2:
    tm_min, tm_max = st.slider("Tm-Bereich (°C)", 48, 75, (58, 62))
    prod_min, prod_max = st.slider("Produktgröße (bp)", 60, 1500, (100, 400))

with c3:
    monoval = st.number_input("Na⁺/K⁺ (mM)", 0.0, 500.0, 50.0, step=1.0)
    dival = st.number_input("Mg²⁺ (mM)", 0.0, 10.0, 1.5, step=0.1)
    dntp = st.number_input("dNTP (mM)", 0.0, 5.0, 0.6, step=0.1)

    st.subheader("🧩 Erweiterte Optionen")
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

    # qPCR Probe
    st.subheader("🧫 qPCR / Probe")
    enable_probe = st.checkbox("Probe mitentwerfen", False)
    reporter = st.selectbox("Reporter", ["FAM","HEX","VIC","ROX","Cy5"], index=0)
    quencher = st.selectbox("Quencher", ["BHQ1","BHQ2","TAMRA","Iowa Black FQ"], index=0)

    # -------------------- Aktionen -------------------------------------------
    run_auto = st.button("🚀 Automatisches Design starten (primer3)")
    st.markdown("---")

    # ================== Automatisches Design =================================
    if run_auto:
        if not P3_OK:
            st.error("❌ `primer3` nicht installiert. Bitte `primer3-py` zur requirements.txt hinzufügen.")
            st.stop()
        if not target_seq and not consensus:
            st.warning("Bitte Zielsequenz eingeben **oder** Alignment für Konsensus laden.")
            st.stop()

        template = consensus if (allow_degenerate and consensus) else target_seq
        if not template:
            st.warning("Keine gültige Vorlage gefunden.")
            st.stop()

        # Extensions bestimmen
        extL, extR = EXT_PRESETS.get(preset, ("",""))
        if custom_left:  extL = custom_left.upper().replace("U","T")
        if custom_right: extR = custom_right.upper().replace("U","T")

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
            {"SEQUENCE_ID": "target", "SEQUENCE_TEMPLATE": template},
            args
        )

        pairs = []
        for i in range(design.get("PRIMER_PAIR_NUM_RETURNED", 0)):
            lseq = design.get(f"PRIMER_LEFT_{i}_SEQUENCE")
            rseq = design.get(f"PRIMER_RIGHT_{i}_SEQUENCE")
            if not lseq or not rseq: 
                continue

            # Optional 5'-Extensions hinzufügen
            lseq_full = (extL + lseq)
            rseq_full = (extR + rseq)

            tm_l = design.get(f"PRIMER_LEFT_{i}_TM")
            tm_r = design.get(f"PRIMER_RIGHT_{i}_TM")
            gc_l, gc_r = gc_percent(lseq_full), gc_percent(rseq_full)

            # Thermodynamik
            dghp = min(dG_hairpin(lseq_full), dG_hairpin(rseq_full))
            dgself = min(dG_homodimer(lseq_full), dG_homodimer(rseq_full))
            dgcross = dG_heterodimer(lseq_full, rseq_full)

            homo = max(longest_homopolymer(lseq_full), longest_homopolymer(rseq_full))
            if homo > max_homopoly:
                continue
            if gc_clamp and (lseq_full[-1] not in "GC" or rseq_full[-1] not in "GC"):
                continue

            # Off-target Screening (optional DB)
            off_hits = 0
            if offtarget_db:
                off_hits = count_offtargets(lseq_full, offtarget_db, mismatches) + \
                           count_offtargets(rseq_full, offtarget_db, mismatches)

            score = simple_score(tm_l, tm_r, gc_l, gc_r, dgself, dgcross, homo, (tm_min+tm_max)/2)

            prod_size = design.get(f"PRIMER_PAIR_{i}_PRODUCT_SIZE")
            lpos = design.get(f"PRIMER_LEFT_{i}")[0]
            rpos = design.get(f"PRIMER_RIGHT_{i}")[0]

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
                "OffTarget_hits": int(off_hits),
                "Amplicon_bp": int(prod_size),
                "Left_start": int(lpos),
                "Right_start": int(rpos),
                "Score": score
            })

        if not pairs:
            st.warning("Keine geeigneten Primerpaare gefunden – Parameter lockern?")
            st.stop()

        df = pd.DataFrame(pairs).sort_values(["OffTarget_hits","Score"], ascending=[True, False]).reset_index(drop=True)
        st.success(f"✅ {len(df)} Primerpaare gefunden")
        st.dataframe(df, use_container_width=True)

        # Amplicon & Probe
        best = df.iloc[0]
        amp_start = best.Left_start + len(best.Left) - len(extL)
        amp_end   = best.Right_start  # inkl. rechter Primerlänge in Plot unten

        # Visualisierung
        if target_seq:
            st.subheader("🧬 Amplicon-Vorschau")
            fig, ax = plt.subplots(figsize=(8,1.8))
            ax.set_xlim(0, len(template))
            ax.set_ylim(0,1)
            ax.set_yticks([]); ax.set_xlabel("Position (bp)")
            ax.axvspan(best.Left_start, best.Left_start+len(best.Left), ymin=0.6, ymax=0.95, color="tab:blue", alpha=0.35, label="Left Primer")
            ax.axvspan(best.Right_start-len(best.Right)+1, best.Right_start+1, ymin=0.6, ymax=0.95, color="tab:orange", alpha=0.35, label="Right Primer")
            ax.axvspan(amp_start, best.Right_start-len(best.Right)+1, ymin=0.25, ymax=0.55, color="tab:green", alpha=0.25, label="Amplicon")
            ax.legend(loc="upper right", ncol=3, fontsize=8, frameon=False)
            st.pyplot(fig, use_container_width=True)

            amplicon_seq = template[amp_start: best.Right_start-len(best.Right)+1]
            st.code(textwrap.fill(amplicon_seq, 80), language="text")
            st.download_button("⬇️ Amplicon als FASTA", f">amplicon\n{amplicon_seq}\n", file_name="amplicon.fasta")

        # qPCR-Probe (einfaches Center-Heuristik-Design)
        if enable_probe and target_seq:
            st.subheader("🧫 qPCR-Probe (heuristisch)")
            center = (amp_start + (best.Right_start-len(best.Right)+1))//2
            p_len = min(28, max(18, int(0.4*(best.Amplicon_bp))))
            p_start = max(0, center - p_len//2)
            probe = template[p_start:p_start+p_len]
            st.write(f"**Probe:** {probe}")
            st.caption(f"Reporter: **{reporter}**, Quencher: **{quencher}**")

        # Export
        st.download_button("⬇️ Primerpaare (CSV)", df.to_csv(index=False).encode("utf-8"),
                           file_name="primer_pairs_advanced.csv", mime="text/csv")

    # ================== Manuelles Prüfen / Scoring ===========================
    st.markdown("---")
    st.subheader("✍️ Manuelles Primer-Scoring")
    m1, m2 = st.columns(2)
    with m1:
        man_left = st.text_input("Left Primer (5'→3')", "")
    with m2:
        man_right = st.text_input("Right Primer (5'→3')", "")

    if st.button("🧮 Manuelle Bewertung"):
        if not man_left or not man_right:
            st.warning("Bitte beide Primer eingeben.")
        else:
            L, R = man_left.upper(), man_right.upper()
            tm_l = primer3.calcTm(L)
            tm_r = primer3.calcTm(R)
            score = simple_score(tm_l, tm_r, gc_percent(L), gc_percent(R),
                                 dG_homodimer(L), dG_heterodimer(L,R),
                                 max(longest_homopolymer(L), longest_homopolymer(R)),
                                 (tm_min+tm_max)/2)
            st.write(f"**Tm**: {tm_l:.1f} / {tm_r:.1f} °C · **GC**: {gc_percent(L):.1f}% / {gc_percent(R):.1f}%")
            st.write(f"**ΔG self/cross (kcal/mol)**: {dG_homodimer(L):.2f} / {dG_heterodimer(L,R):.2f}")
            st.success(f"Score (0..1): **{score}**")

    # ================== AI-Zusammenfassung (optional) ========================
    st.markdown("---")
    st.subheader("🤖 AI-Zusammenfassung (optional)")
    ai_btn = st.button("Ergebnisse erklären lassen")
    if ai_btn:
        summary = (
            "Dieses Modul entwirft Primerpaare mittels primer3, prüft Hairpins/Dimer (ΔG), "
            "berechnet GC/Tm, erlaubt 5'-Extensions, degenerate Konsensus-Primer aus Alignments, "
            "führt eine einfache Off-Target-Suche gegen lokale FASTA durch und extrahiert das Amplicon. "
            "qPCR-Proben können heuristisch erstellt werden."
        )
        if OPENAI_OK and os.getenv("OPENAI_API_KEY"):
            try:
                openai.api_key = os.getenv("OPENAI_API_KEY")
                prompt = ("Erkläre einem Laboranwender knapp auf Deutsch, wie die erzeugten Primer bewertet "
                          "werden (Tm, GC, ΔG, Off-Target-Hits, Amplicon-Größe) und gib Tipps zur Optimierung.")
                resp = openai.ChatCompletion.create(
                    model="gpt-4o-mini",
                    messages=[{"role":"system","content":prompt}]
                )
                summary = resp["choices"][0]["message"]["content"]
            except Exception:
                pass
        st.info(summary)
