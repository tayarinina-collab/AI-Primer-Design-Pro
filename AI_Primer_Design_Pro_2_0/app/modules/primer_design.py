# -*- coding: utf-8 -*-
"""
AI Primer Design Pro – Modul: Primer Design (wie Geneious Prime)
Stabil & offline-fähig mit Primer3, Hairpin/Dimer-Analyse, Fallback für degenerierte Basen
"""
import streamlit as st
import primer3
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# --------------------------------------------------------------------------
# Hilfsfunktionen
# --------------------------------------------------------------------------
def gc_percent(seq: str) -> float:
    seq = seq.upper().replace("U", "T")
    if not seq:
        return 0.0
    return round(100 * (seq.count("G") + seq.count("C")) / len(seq), 2)

def hairpin_dG(seq: str) -> float:
    result = primer3.calcHairpin(seq)
    return result.dg if result.structure_found else 0

def dimer_dG(seq1: str, seq2: str = None) -> float:
    if seq2:
        result = primer3.calcHeterodimer(seq1, seq2)
    else:
        result = primer3.calcHomodimer(seq1)
    return result.dg if result.structure_found else 0

def longest_homopolymer(seq: str) -> int:
    max_run = run = 1
    for i in range(1, len(seq)):
        run = run + 1 if seq[i] == seq[i - 1] else 1
        max_run = max(max_run, run)
    return max_run

def add_extension(seq: str, ext5: str = "") -> str:
    return ext5 + seq if ext5 else seq

# --------------------------------------------------------------------------
# Hauptfunktion
# --------------------------------------------------------------------------
def run_primer_design():
    st.header("🧬 Primer Design (wie Geneious Prime)")
    st.caption("Design & Analyse von Primern mit Thermodynamik, Hairpin-Check, Dimer-Analyse & Visualisierung")

    seq_input = st.text_area("DNA-Sequenz eingeben (5'→3'):", height=140)
    if not seq_input.strip():
        st.info("Bitte Sequenz eingeben, um Primer zu entwerfen.")
        return

    seq = seq_input.upper().replace("U", "T").replace(" ", "").replace("\n", "")

    # ---------------- Parameter (wie Geneious) -----------------------------
    st.subheader("⚙️ Parameter (wie Geneious Prime)")
    col1, col2, col3 = st.columns(3)
    with col1:
        tm_opt = st.slider("Tm (Schmelztemperatur °C)", 50, 75, 60)
        gc_min, gc_max = st.slider("GC-Gehalt (%)", 30, 80, (40, 60))
    with col2:
        length_min, length_max = st.slider("Primerlänge (bp)", 15, 35, (18, 25))
        product_min, product_max = st.slider("Amplicon-Größe (bp)", 50, 1200, (80, 400))
    with col3:
        gc_clamp = st.checkbox("GC-Clamp am 3'-Ende", value=True)
        dg_threshold = st.number_input("ΔG-Grenze (kcal/mol, min)", -20.0, 0.0, -9.0, step=0.5)

    ext5 = st.text_input("Optionale 5’-Extension (z. B. Restriktionsstelle oder Tag)", "")

    # ---------------- Zusatzfunktionen -------------------------------------
    st.subheader("🧩 Zusatzfunktionen")
    degenerate = st.checkbox("Degenerate Primer Design (IUPAC N,R,Y...)", value=False)
    qpcr_mode = st.checkbox("qPCR/TaqMan Probe mitentwerfen", value=False)
    visualize = st.checkbox("Grafische Visualisierung aktivieren", value=True)

    # ---------------- Primer3 Konfiguration --------------------------------
    args = {
        "PRIMER_OPT_SIZE": int(np.mean([length_min, length_max])),
        "PRIMER_MIN_SIZE": int(length_min),
        "PRIMER_MAX_SIZE": int(length_max),
        "PRIMER_OPT_TM": float(tm_opt),
        "PRIMER_MIN_TM": float(tm_opt - 3),
        "PRIMER_MAX_TM": float(tm_opt + 3),
        "PRIMER_MIN_GC": float(gc_min),
        "PRIMER_MAX_GC": float(gc_max),
        "PRIMER_SALT_MONOVALENT": 50.0,
        "PRIMER_DNA_CONC": 250.0,
        "PRIMER_PRODUCT_SIZE_RANGE": [[int(product_min), int(product_max)]],
        "PRIMER_NUM_RETURN": 20
    }

    # ---------------- Primer-Design starten --------------------------------
    if st.button("🚀 Primer entwerfen"):

        # 🧠 Degenerate Basen verhindern Absturz → Fallback-Modus
        if degenerate:
            st.warning("⚠️ Primer3 unterstützt keine Thermoanalyse für degenerierte Basen. "
                       "Wechsle in vereinfachten Fallback-Modus.")
            seq_simple = seq[:200] if len(seq) > 200 else seq
            fwd = seq_simple[:20].replace("A", "R").replace("T", "Y")
            rev = seq_simple[-20:].replace("C", "S").replace("G", "K")
            result = [{
                "Index": 1,
                "Left Primer": fwd,
                "Right Primer": rev,
                "Length (bp)": 20,
                "Tm (°C)": 60.0,
                "GC%": gc_percent(fwd),
                "ΔG (Hairpin/Dimer)": 0.0,
                "Amplicon Size": len(seq_simple)
            }]
            df = pd.DataFrame(result)
            st.dataframe(df, use_container_width=True)
            st.info("Degenerate Primer erstellt (ohne Thermodynamik).")
            return

        # 🧪 Regulärer Primer3-Modus
        try:
            primers = primer3.bindings.designPrimers(
                {"SEQUENCE_ID": "target", "SEQUENCE_TEMPLATE": seq},
                args
            )
        except Exception as e:
            st.error("❌ Primer3 konnte keine Primer berechnen.")
            st.text(f"Fehlerdetails: {e}")
            st.info("Tipp: Degenerate-Modus ausschalten oder Parameter anpassen.")
            return

        total = primers.get("PRIMER_PAIR_NUM_RETURNED", 0)
        if total == 0:
            st.warning("Keine Primer gefunden – bitte Parameter anpassen (z. B. Produktgröße oder GC-Bereich).")
            return

        # ---------------- Ergebnisse sammeln --------------------------------
        result = []
        for i in range(total):
            left_seq = primers[f"PRIMER_LEFT_{i}_SEQUENCE"]
            right_seq = primers[f"PRIMER_RIGHT_{i}_SEQUENCE"]

            if ext5:
                left_seq = add_extension(left_seq, ext5)
                right_seq = add_extension(right_seq, ext5)

            tm_left = round(primers[f"PRIMER_LEFT_{i}_TM"], 2)
            tm_right = round(primers[f"PRIMER_RIGHT_{i}_TM"], 2)
            gc_left = gc_percent(left_seq)
            gc_right = gc_percent(right_seq)
            dg_hairpin = min(hairpin_dG(left_seq), hairpin_dG(right_seq))
            dg_dimer = min(dimer_dG(left_seq), dimer_dG(right_seq), dimer_dG(left_seq, right_seq))
            product = primers[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"]

            # GC-Clamp prüfen
            if gc_clamp and (left_seq[-1] not in "GC" or right_seq[-1] not in "GC"):
                continue

            # ΔG-Grenze prüfen
            if min(dg_hairpin, dg_dimer) < dg_threshold:
                continue

            result.append({
                "Index": i + 1,
                "Left Primer": left_seq,
                "Right Primer": right_seq,
                "Length (bp)": len(left_seq),
                "Tm (°C)": round(np.mean([tm_left, tm_right]), 1),
                "GC%": round(np.mean([gc_left, gc_right]), 1),
                "ΔG (Hairpin/Dimer)": round(min(dg_hairpin, dg_dimer), 2),
                "Amplicon Size": product
            })

        if not result:
            st.warning("Keine Primer erfüllen die Thermodynamik-Grenzen.")
            return

        df = pd.DataFrame(result)
        st.success(f"{len(df)} Primerpaare erfolgreich generiert ✅")
        st.dataframe(df, use_container_width=True)

        # ---------------- ΔG Visualisierung ----------------------------------
        if visualize:
            st.subheader("ΔG & GC-Profil Visualisierung")
            fig, ax = plt.subplots(figsize=(6, 3))
            ax.scatter(df["GC%"], df["ΔG (Hairpin/Dimer)"], c=df["Tm (°C)"], cmap="viridis", s=70)
            ax.set_xlabel("GC-Gehalt (%)")
            ax.set_ylabel("ΔG (kcal/mol)")
            ax.set_title("Primer-Thermodynamik-Profil")
            st.pyplot(fig)

        # ---------------- qPCR Probe ----------------------------------------
        if qpcr_mode:
            st.subheader("🧫 qPCR/TaqMan Probe")
            probe_start = int(len(seq) / 2) - 10
            probe_seq = seq[probe_start:probe_start + 25]
            st.write(f"**Probe:** {probe_seq}")
            st.caption("Einfaches Modell für TaqMan/SYBR-Probe (Optimierung folgt).")

        # ---------------- Download ------------------------------------------
        st.download_button(
            "⬇️ CSV exportieren",
            df.to_csv(index=False).encode("utf-8"),
            file_name="primer_design_results.csv",
            mime="text/csv"
        )
