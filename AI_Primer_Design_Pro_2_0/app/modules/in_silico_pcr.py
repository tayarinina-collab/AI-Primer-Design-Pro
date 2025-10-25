import streamlit as st
from Bio.Seq import Seq

def simple_in_silico_pcr(template: str, fwd: str, rev: str, max_product=3000):
    template = template.upper()
    fwd = fwd.upper()
    rev_rc = str(Seq(rev).reverse_complement()).upper()
    start = template.find(fwd)
    end = template.find(rev_rc)
    if start != -1 and end != -1 and end > start:
        product = template[start:end+len(rev_rc)]
        if len(product) <= max_product:
            return product
    return None

def render():
    st.header("In-Silico PCR")
    st.caption("Virtual amplification with simple binding search / Virtuelle Amplifikation mit einfacher Bindungssuche")

    template = st.text_area("Template sequence (5'->3')", height=160, value=("ATGCGTACGTAGCTAGCTAGCTAGCTAGCATCGATCGATGCTAGCTAGCTAGC"*3))
    c1, c2 = st.columns(2)
    with c1:
        fwd = st.text_input("Forward primer (5'->3')", value="ATGCGTACGTAGC")
    with c2:
        rev = st.text_input("Reverse primer (5'->3')", value="GCTAGCTAGCATC")

    max_product = st.number_input("Max product size (bp)", min_value=100, value=2000, step=50)

    if st.button("Run In-Silico PCR / PCR simulieren", use_container_width=True):
        prod = simple_in_silico_pcr(template, fwd, rev, max_product=max_product)
        if prod:
            st.success(f"Amplified product: {len(prod)} bp")
            st.code(prod[:500] + ("..." if len(prod)>500 else ""))
        else:
            st.warning("No valid amplification found with given primers.")