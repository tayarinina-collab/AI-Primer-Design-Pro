import streamlit as st

def render():
    st.header("Einstellungen / Info")
    st.markdown("""
AI Primer Design Pro v2.0 (Beta)

- Bilingual UI (DE/EN)
- Core modules: Primer Design, In-Silico PCR, Plasmid Map, Restriction Tools, Protein Tools

Roadmap:
- Databases (NCBI/UniProt), 3D Viewer, MSA/Phylogeny, Reports (PDF), Workflow Automation

Tips:
- Use shorter demo sequences for speed in the cloud
- Add OPENAI_API_KEY in Streamlit Secrets for the chatbot

Credits:
Built with Streamlit, BioPython, primer3-py, Plotly.
""", unsafe_allow_html=False)