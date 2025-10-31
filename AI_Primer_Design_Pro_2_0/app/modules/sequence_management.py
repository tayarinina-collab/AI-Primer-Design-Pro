import streamlit as st
from Bio.SeqUtils import gc_fraction
from Bio.Seq import Seq
import pandas as pd
import matplotlib.pyplot as plt
from io import StringIO

def detect_sequence_type(seq: str) -> str:
    dna_bases = set("ATGCN")
    rna_bases = set("AUGCN")
    protein_letters = set("ACDEFGHIKLMNPQRSTVWY")

    seq_clean = "".join(seq.upper().split())  # entfernt Leerzeichen & Zeilenumbrüche
    seq_set = set(seq_clean)

    if seq_set <= dna_bases:
        return "DNA"
    elif seq_set <= rna_bases:
        return "RNA"
    elif seq_set <= protein_letters:
        return "Protein"
    else:
        return "Unbekannt"

def run_sequence_management():
    st.header("🧬 Sequence Management")
    st.caption("Analyse und Verwaltung von DNA-, RNA- oder Protein-Sequenzen")

    uploaded_file = st.file_uploader("Datei hochladen (FASTA, GenBank, TXT)", type=["fasta", "txt", "gb"])
    text_input = st.text_area("Oder Sequenz direkt einfügen:", height=150)

    sequence_data = ""
    if uploaded_file:
        sequence_data = uploaded_file.getvalue().decode("utf-8")
    elif text_input.strip():
        sequence_data = text_input.strip()

    if sequence_data:
        st.subheader("🔍 Sequenzanalyse")
        seq_str = sequence_data.replace("\n", "").replace(" ", "").upper()
        seq = Seq(seq_str)

        # 🧩 präzise Typ-Erkennung
        seq_type = detect_sequence_type(seq_str)
        length = len(seq)
        gc_content = round(gc_fraction(seq) * 100, 2) if seq_type in ["DNA", "RNA"] else "N/A"

        # Anzeige der Basisinformationen
        st.write(f"**Sequenztyp:** {seq_type}")
        st.write(f"**Länge:** {length} bp / aa")
        st.write(f"**GC-Gehalt:** {gc_content}%")

        # GC-Profilplot (nur DNA/RNA)
        if seq_type in ["DNA", "RNA"] and length >= 50:
            window = 50
            gc_profile = [round(gc_fraction(seq[i:i+window]) * 100, 2)
                          for i in range(0, len(seq)-window, window)]
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
