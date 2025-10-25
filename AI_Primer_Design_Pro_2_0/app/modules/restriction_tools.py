import streamlit as st
from Bio.Seq import Seq
from Bio.Restriction import RestrictionBatch, Analysis, AllEnzymes

def render():
    st.header("Restriktions-Analyse")
    st.caption("Find cut sites using Bio.Restriction / Schnittstellen mit Bio.Restriction finden")

    seq = st.text_area("DNA sequence (5'->3')", height=150, value="GAATTCGCGGCCGCAAGCTT" * 5)
    preset = st.selectbox("Enzyme preset", ["Common enzymes", "All enzymes"], index=0)
    if preset == "Common enzymes":
        enzymes = ["EcoRI", "NotI", "HindIII", "BamHI", "XhoI", "NheI", "XbaI"]
    else:
        enzymes = [e.__name__ for e in AllEnzymes]

    custom = st.text_input("Custom enzyme list (comma separated)", value="")
    if custom.strip():
        enzymes = [e.strip() for e in custom.split(",") if e.strip()]

    if st.button("Analyze / Analysieren", use_container_width=True):
        rb = RestrictionBatch(enzymes)
        ana = Analysis(rb, Seq(seq.upper()))
        result = ana.full()
        if not result:
            st.warning("No cut sites found for selected enzymes.")
        else:
            for enz, sites in result.items():
                st.write(f"{enz}: {sites}")