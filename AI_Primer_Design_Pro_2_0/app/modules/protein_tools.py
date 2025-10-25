import streamlit as st
from Bio.SeqUtils.ProtParam import ProteinAnalysis

def render():
    st.header("Protein-Tools + KI-Assistent")
    st.caption("Physicochemical properties + optional AI helper / Physikochemische Eigenschaften + optionaler KI-Assistent")

    seq = st.text_area("Protein sequence (1-letter codes)", height=140, value="MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQR")
    if st.button("Analyze / Analysieren"):
        try:
            pa = ProteinAnalysis(seq.replace("\n","").upper())
            st.write({
                "Length": len(seq),
                "MW": round(pa.molecular_weight(), 2),
                "pI": round(pa.isoelectric_point(), 2),
                "Aromaticity": round(pa.aromaticity(), 3)
            })
            st.line_chart(pa.protein_scale(window=7, param_dict=ProteinAnalysis.protein_scale_params['kd']))
        except Exception as e:
            st.error(str(e))

    st.subheader("AI Assistant (optional)")
    st.info("Set environment variable OPENAI_API_KEY in Streamlit Cloud Secrets for the chatbot. This demo will return a template if no key is set.")

    question = st.text_input("Ask a protein design question / Frage stellen")
    if st.button("Ask / Fragen"):
        import os
        api_key = os.environ.get("OPENAI_API_KEY")
        if not api_key:
            st.warning("No OPENAI_API_KEY detected. Returning a template answer.")
            st.write("Template answer: Consider pI, hydrophobicity, secondary structure, and avoid primer dimers/hairpins.")
        else:
            st.success("OPENAI_API_KEY detected. Add your API call in this function.")