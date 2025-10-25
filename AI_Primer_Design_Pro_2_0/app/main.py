import streamlit as st
from modules.sequence_management import run_sequence_management
from modules.primer_design import run_primer_design
from modules.in_silico_pcr import run_in_silico_pcr
from modules.protein_tools import run_protein_tools
from modules.plasmid_designer import run_plasmid_designer
from modules.ui_layout import set_theme

# --- Sidebar Navigation ---
st.set_page_config(page_title="AI Primer Design Pro", layout="wide", page_icon="🧬")

set_theme()  # Dark/Light Theme Setup

st.sidebar.title("🧬 AI Primer Design Pro")
st.sidebar.markdown("**Intelligente Bioinformatik-Plattform für moderne Labore**")

# Navigation: Module-Übersicht
menu = st.sidebar.radio(
    "🧩 Module auswählen / Select Module",
    [
        "🏠 Übersicht",
        "🧬 Sequence Management",
        "🧫 Primer Design",
        "🧪 In-Silico PCR",
        "🧫 Protein Tools",
        "🧫 Plasmid Designer"
    ],
)

# --- Hauptansicht ---
if menu == "🏠 Übersicht":
    st.title("Willkommen in AI Primer Design Pro 🧬")
    st.markdown("""
    **Deutsch 🇩🇪**  
    Willkommen bei *AI Primer Design Pro*!  
    Diese Plattform kombiniert **Bioinformatik**, **KI-Analyse** und **modernes Labor-Design**,  
    um deine molekularbiologischen Workflows zu automatisieren.

    **English 🇬🇧**  
    Welcome to *AI Primer Design Pro*!  
    This platform unites **bioinformatics**, **AI analysis**, and **modern lab design**  
    for seamless automation of molecular biology workflows.
    """)

# --- Modulaufrufe ---
elif menu == "🧬 Sequence Management":
    run_sequence_management()

elif menu == "🧫 Primer Design":
    run_primer_design()

elif menu == "🧪 In-Silico PCR":
    run_in_silico_pcr()

elif menu == "🧫 Protein Tools":
    run_protein_tools()

elif menu == "🧫 Plasmid Designer":
    run_plasmid_designer()
