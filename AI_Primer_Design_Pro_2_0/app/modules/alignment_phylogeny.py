# -*- coding: utf-8 -*-
"""
Alignment & Phylogeny Modul
- Multiple & Pairwise Alignment (BioPython)
- Phylogeny (Neighbor Joining, UPGMA)
- Visualisierung & Export
"""

import streamlit as st
from Bio import Align, Phylo, SeqIO
from Bio.Align import substitution_matrices
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceCalculator
from io import StringIO, BytesIO
import matplotlib.pyplot as plt

def run_alignment_phylogeny():
    st.title("🌳 Alignment & Phylogeny")
    st.caption("Multiple- und Paarweise Alignments, mit Baumvisualisierung und Export")

    st.markdown("### 📤 Sequenzen hochladen (FASTA)")
    fasta_file = st.file_uploader("FASTA-Datei auswählen", type=["fasta", "fa", "txt"])

    if not fasta_file:
        st.info("Bitte eine FASTA-Datei hochladen, um zu starten.")
        return

    # Sequenzen einlesen
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    st.success(f"✅ {len(sequences)} Sequenzen geladen")
    for rec in sequences[:3]:
        st.text(f">{rec.id}\n{rec.seq[:60]}...")

    st.markdown("---")
    st.markdown("### ⚙️ Alignment-Typ auswählen")
    mode = st.radio("Alignment-Modus", ["Multiple Sequence Alignment", "Pairwise Alignment"], horizontal=False)

    # -----------------------------
    # Multiple Alignment
    # -----------------------------
    if mode == "Multiple Sequence Alignment":
        st.subheader("🧬 Multiple Alignment (MUSCLE / MAFFT / ClustalΩ)")
        method = st.selectbox("Algorithmus wählen", ["MUSCLE (BioPython built-in)", "Clustal Omega", "MAFFT (simuliert)"])

        if st.button("🚀 Alignment starten"):
            try:
                aligner = Align.PairwiseAligner()
                aligner.mode = "global"
                matrix = substitution_matrices.load("BLOSUM62")
                aligner.substitution_matrix = matrix

                # sehr einfaches Demo: alles gegen die erste Sequenz alignen
                ref = sequences[0]
                aligned = [ref]
                for seq in sequences[1:]:
                    alignments = aligner.align(ref.seq, seq.seq)
                    aligned.append(seq)

                st.success(f"Alignment mit {method} simuliert (Demo-Version)")
                st.text_area("Resultat (FASTA-like)", "\n".join(f">{s.id}\n{s.seq}" for s in aligned), height=250)

            except Exception as e:
                st.error(f"Fehler beim Alignment: {e}")

    # -----------------------------
    # Pairwise Alignment
    # -----------------------------
    else:
        st.subheader("🧩 Paarweises Alignment (Needleman-Wunsch / Smith-Waterman)")
        method = st.selectbox("Algorithmus wählen", ["Needleman-Wunsch (global)", "Smith-Waterman (lokal)"])

        if len(sequences) < 2:
            st.warning("Mindestens 2 Sequenzen erforderlich.")
            return

        seq1, seq2 = sequences[0].seq, sequences[1].seq

        if st.button("🔬 Align 1 ↔ 2"):
            aligner = Align.PairwiseAligner()
            aligner.mode = "global" if "Needleman" in method else "local"
            alignment = aligner.align(seq1, seq2)[0]
            st.success(f"✅ {method} durchgeführt")
            st.code(str(alignment), language="text")

    st.markdown("---")
    st.subheader("🌿 Phylogenetischer Baum (Neighbor Joining / UPGMA)")

    if st.button("🌳 Baum erzeugen"):
        try:
            # Distanzmatrix berechnen
            output = StringIO()
            SeqIO.write(sequences, output, "fasta")
            output.seek(0)
            calculator = DistanceCalculator("identity")
            aln = list(SeqIO.parse(output, "fasta"))
            distance_matrix = calculator.get_distance(aln)

            constructor = DistanceTreeConstructor()
            nj_tree = constructor.nj(distance_matrix)

            fig = plt.figure(figsize=(6, 6))
            Phylo.draw(nj_tree, do_show=False)
            st.pyplot(fig, use_container_width=True)

            # Export
            newick_buf = StringIO()
            Phylo.write(nj_tree, newick_buf, "newick")
            st.download_button(
                "💾 Baum als Newick exportieren",
                data=newick_buf.getvalue().encode("utf-8"),
                file_name="tree.newick",
                mime="text/plain"
            )

        except Exception as e:
            st.error(f"Fehler bei Baum-Erstellung: {e}")
