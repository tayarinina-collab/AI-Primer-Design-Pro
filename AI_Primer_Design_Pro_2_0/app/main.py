from modules.protein_tools import run_protein_tools
from modules.cloning_tools import run_cloning_tools
from modules.primer_design import run_primer_design
from modules.primer_design_advanced import run_primer_design_advanced
import streamlit as st
# --- Importiere alle Module ---
from modules.sequence_management import run_sequence_management
from modules.primer_design import run_primer_design
from modules.in_silico_pcr import run_in_silico_pcr
from modules.protein_tools import run_protein_tools
from modules.plasmid_designer import run_plasmid_designer
from modules.ui_layout import set_theme

# --- Seiteneinstellungen ---
st.set_page_config(
    page_title="AI Primer Design Pro",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# --- Theme Setup ---
set_theme()  # Schaltet zwischen Light und Dark Mode

# --- Sidebar Titel ---
st.sidebar.title("🧬 AI Primer Design Pro")
st.sidebar.markdown("**Intelligente Bioinformatik-Plattform für moderne Labore**")
st.sidebar.markdown("---")

# --- Sprachumschalter ---
language = st.sidebar.radio("🌍 Sprache / Language", ["🇩🇪 Deutsch", "🇬🇧 English"])

st.sidebar.markdown("---")

# --- Navigationsmenü ---
menu = st.sidebar.radio(
    "🧩 Module auswählen / Select Module",
    [
        "🏠 Übersicht",
        "🧬 Sequence Management",
        "🧫 Primer Design",
        "🧪 Primer Design – Advanced",
        "🧫 Cloning & Assembly Tools",
        "🧬 Protein Tools",
        "🧫 Plasmid Designer"
    ],
)

# --- HAUPTINHALT ---
if menu == "🏠 Übersicht":
    if language == "🇩🇪 Deutsch":
        st.title("Willkommen in AI Primer Design Pro 🧬")
        st.markdown("""
        Willkommen bei **AI Primer Design Pro**,  
        deiner intelligenten Bioinformatik-Plattform für DNA-, RNA- und Protein-Analysen.  
        Hier kombinieren sich **KI**, **Laborautomatisierung** und **visuelle Werkzeuge**,  
        um Forschungsprozesse zu vereinfachen und zu beschleunigen.  
        """)
        st.info("🌗 Tipp: Du kannst im Seitenmenü zwischen **Dark- und Light-Mode** wechseln.")
    else:
        st.title("Welcome to AI Primer Design Pro 🧬")
        st.markdown("""
        Welcome to **AI Primer Design Pro**,  
        your intelligent bioinformatics platform for DNA, RNA, and protein analysis.  
        Combining **AI**, **automation**, and **visual lab tools**  
        to simplify and accelerate research workflows.  
        """)
        st.info("🌗 Tip: You can switch between **Dark and Light mode** in the sidebar.")

# --- SEQUENCE MANAGEMENT ---
elif menu == "🧬 Sequence Management":
    run_sequence_management()

# --- PRIMER DESIGN ---
elif menu == "🧫 Primer Design":
    run_primer_design()
elif menu == "🧪 Primer Design – Advanced":
    run_primer_design_advanced()

# --- IN-SILICO PCR ---
elif menu == "🧫 Cloning & Assembly Tools":
    run_cloning_tools()

# --- PROTEIN TOOLS ---
elif menu == "🔬 Protein Tools":
    run_protein_tools()

# --- PLASMID DESIGNER ---
elif menu == "🧫 Plasmid Designer":
    run_plasmid_designer()

# --- Footer ---
st.markdown("---")
st.caption("🧠 Entwickelt mit ❤️ in Hamburg · Version 2.0 · Zweisprachig DE/EN")
