import streamlit as st
from Bio import pairwise2

def render():
    st.header("Alignments (Basic)")
    st.caption("Pairwise alignment demo / Einfacher Pairwise-Alignment-Modus")

    s1 = st.text_area("Sequence 1", value="GATTACA")
    s2 = st.text_area("Sequence 2", value="GCATGCU")

    if st.button("Align"):
        alns = pairwise2.align.globalxx(s1, s2)
        if not alns:
            st.warning("No alignment found.")
        else:
            aln = alns[0]
            st.code(f"{aln.seqA}\n{aln.seqB}\nScore: {aln.score}")
            st.caption("MSA/MAFFT integration planned.")