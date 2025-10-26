# -*- coding: utf-8 -*-
"""
Cloning & Assembly Tools (Geneious Pro Style)
Offline-Module: Gibson, Golden Gate, Restriktionsanalyse, Plasmidmap
"""
import streamlit as st
import pandas as pd
from Bio.Seq import Seq
from Bio import Restriction
from Bio.Restriction import *
import io
import plotly.graph_objects as go
from Bio import SeqIO

# ============================ UI ===================================
def run_cloning_tools():
    st.title("🧫 Cloning & Assembly Tools")
    st.caption("Gibson, Golden Gate, Restriktionsenzyme, In-Silico Verdau & Ligation, Plasmid Maps")

    tabs = st.tabs(["🔬 Gibson Assembly", "🧩 Golden Gate", "✂️ Restriktionsanalyse", "🧬 In-Silico Verdau", "🧭 Plasmid Viewer"])

    # ---------- Gibson Assembly ----------
    with tabs[0]:
        st.subheader("🔬 Gibson Assembly Assistant")
        seq1 = st.text_area("Fragment 1 (5'→3')", height=120)
        seq2 = st.text_area("Fragment 2 (5'→3')", height=120)
        overlap = st.slider("Overlap (bp)", 10, 60, 20)

        if st.button("🧠 Gibson Assembly simulieren"):
            if seq1 and seq2:
                overlap_seq = seq1[-overlap:]
                if overlap_seq in seq2[:overlap]:
                    assembled = seq1 + seq2[overlap:]
                    st.success("✅ Gibson erfolgreich! Overlap erkannt.")
                    st.code(assembled)
                else:
                    st.warning("⚠️ Kein passender Overlap gefunden – überprüfe Sequenzen.")
            else:
                st.info("Bitte beide Sequenzen eingeben.")

    # ---------- Golden Gate ----------
    with tabs[1]:
        st.subheader("🧩 Golden Gate Assembly Designer")
        seq = st.text_area("DNA-Sequenz für Golden Gate", height=120)
        enzyme = st.selectbox("Typ IIS Restriktionsenzym", ["BsaI", "BsmBI", "BbsI", "Esp3I"])

        if st.button("🔧 Simulation starten"):
            st.info(f"Simuliere Golden Gate mit {enzyme} ... (Schnittstellen werden gesucht)")
            # einfache Erkennung von Typ IIS Motiven:
            motif = {"BsaI": "GGTCTC", "BsmBI": "CGTCTC", "BbsI": "GAAGAC", "Esp3I": "CGTCTC"}[enzyme]
            pos = [i for i in range(len(seq)) if seq[i:i+6] == motif]
            if pos:
                st.success(f"🔬 {len(pos)} Schnittstellen gefunden bei Position(en): {pos}")
            else:
                st.warning("Keine Schnittstellen gefunden.")

    # ---------- Restriktionsanalyse ----------
    with tabs[2]:
        st.subheader("✂️ Restriktionsenzyme-Analyse")
        uploaded_seq = st.file_uploader("FASTA-Datei hochladen", type=["fasta","fa","txt"])
        if uploaded_seq:
            recs = list(SeqIO.parse(io.StringIO(uploaded_seq.getvalue().decode()), "fasta"))
            seq = recs[0].seq
            analysis = Restriction.Analysis(Restriction.AllEnzymes, seq)
            result = analysis.full()
            st.write(f"Gefundene Schnittstellen: {len(result)} Enzyme erkannt.")
            st.json(result)

    # ---------- In-Silico Verdau ----------
    with tabs[3]:
        st.subheader("🧬 In-Silico Verdau & Ligation")
        seq = st.text_area("DNA-Sequenz für Verdau", height=120)
        enzyme = st.selectbox("Enzym auswählen", ["EcoRI", "BamHI", "HindIII", "NotI"])
        if st.button("🔪 Verdau simulieren"):
            try:
                cut = getattr(Restriction, enzyme)
                cuts = cut.search(Seq(seq))
                if cuts:
                    st.success(f"{enzyme} schneidet an Position(en): {cuts}")
                else:
                    st.warning("Keine Schnittstelle gefunden.")
            except Exception as e:
                st.error(f"Fehler: {e}")

    # ---------- Plasmid Viewer ----------
    with tabs[4]:
        st.subheader("🧭 Plasmid Map Visualisierung")
        seq = st.text_area("DNA-Sequenz oder FASTA", height=120)
        if st.button("🧫 Plasmid anzeigen"):
            if len(seq) > 0:
                length = len(seq)
                fig = go.Figure()
                fig.add_trace(go.Scatterpolar(
                    r=[1]*length,
                    theta=list(range(length)),
                    mode='lines',
                    line=dict(color='lime', width=3),
                    name="DNA Circle"
                ))
                fig.update_layout(showlegend=False, polar=dict(radialaxis=dict(visible=False)))
                st.plotly_chart(fig, use_container_width=True)
            else:
                st.info("Bitte DNA-Sequenz eingeben.")
