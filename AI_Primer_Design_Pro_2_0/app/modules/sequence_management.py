import streamlit as st
from Bio.SeqUtils import GC
from Bio.Seq import Seq
import pandas as pd
import matplotlib.pyplot as plt
from io import StringIO

def run_sequence_management():
    st.header("🧬 Sequence Management")
    st.caption("Analyse und Verwaltung von DNA-, RNA- oder Protein-Sequenzen")

    # Upload oder Texteingabe
    uploaded_file = st.file_uploader("Datei hochladen (FASTA, GenBank, TXT)", type=["fasta", "txt", "gb"])
    text_input = st.text_area("Oder Sequenz direkt einfügen:", height=150)

    sequence_data = ""
    if uploaded_file:
        sequence_data = uploaded_file.getvalue().decode("utf-8")
    elif text_input.strip():
        sequence_data = text_input.strip()

    if sequence_data:
        st.subheader("🔍 Sequenzanalyse")
        seq = Seq(sequence_data.replace("\n", "").replace(" ", "").upper())

        # Basiseigenschaften
        seq_type = "Protein" if set(seq) <= set("ACDEFGHIKLMNPQRSTVWY") else "RNA" if "U" in seq else "DNA"
        length = len(seq)
        gc_content = GC(seq) if seq_type in ["DNA", "RNA"] else "N/A"

        # Anzeige der Basisinformationen
        st.write(f"**Sequenztyp:** {seq_type}")
        st.write(f"**Länge:** {length} bp / aa")
        st.write(f"**GC-Gehalt:** {gc_content}%")

        # GC-Profilplot
        if seq_type in ["DNA", "RNA"]:
            gc_profile = [GC(seq[i:i+50]) for i in range(0, len(seq)-49, 50)]
            fig, ax = plt.subplots()
            ax.plot(gc_profile)
            ax.set_xlabel("Position (bp)")
            ax.set_ylabel("GC (%)")
            ax.set_title("GC-Profil")
            st.pyplot(fig)

        # Exportoptionen
        st.download_button("📥 Exportiere als FASTA", data=str(seq), file_name="sequence.fasta")
        st.success("✅ Analyse abgeschlossen.")
    else:
        st.info("Bitte Sequenz hochladen oder einfügen, um zu starten.")
