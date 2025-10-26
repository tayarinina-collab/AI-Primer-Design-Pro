# -*- coding: utf-8 -*-
import streamlit as st

# --- Module Imports ---
from modules.alignment_phylogeny import run_alignment_phylogeny
from modules.sequence_management import run_sequence_management
from modules.primer_design import run_primer_design
from modules.primer_design_advanced import run_primer_design_advanced
from modules.cloning_tools import run_cloning_tools
from modules.protein_tools import run_protein_tools
from modules.plasmid_designer import run_plasmid_designer
from modules.plasmid_plus import run_plasmid_plus
from modules.database_integration import run_database_integration
from modules.data_management import run_data_management
from modules.ui_layout import set_theme

# --- Seiteneinstellungen ---
st.set_page_config(
    page_title="AI Primer Design Pro",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# --- Theme Setup ---
set_theme()

# --- Sidebar Titel ---
st.sidebar.markdown("<h2 style='text-align:center;'>🧬 AI Primer Design Pro</h2>", unsafe_allow_html=True)
st.sidebar.caption("Intelligente Bioinformatik-Plattform für moderne Labore")
st.sidebar.markdown("---")

# --- Sprachumschalter ---
language = st.sidebar.radio("🌍 Sprache / Language", ["Deutsch", "English"], horizontal=True)
st.sidebar.markdown("---")

# --- Navigation mit visuellen Icons ---
st.sidebar.markdown("### 🧩 Module auswählen / Select Module")
menu = st.sidebar.radio(
    "Modul-Liste",
    [
        "🏠 Übersicht",
        "🧬 Sequence Management",
        "🧫 Primer Design",
        "🧪 Primer Design – Advanced",
        "🧫 Cloning & Assembly Tools",
        "🧬 Protein Tools",
        "🧫 Plasmid Designer",
        "🧬 Plasmid Plus",
        "🧫 Database & Reference Integration",
        "🧬 Data Management"
        "🌳 Alignment & Phylogeny",
    ],
)

# --- HAUPTINHALT ---
if menu == "🏠 Übersicht":
    if language == "Deutsch":
        st.title("Willkommen bei AI Primer Design Pro 🧬")
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

# --- MODULE: Sequence Management ---
elif menu == "🧬 Sequence Management":
    run_sequence_management()

# --- MODULE: Primer Design ---
elif menu == "🧫 Primer Design":
    run_primer_design()

# --- MODULE: Primer Design – Advanced ---
elif menu == "🧪 Primer Design – Advanced":
    run_primer_design_advanced()

# --- MODULE: Cloning & Assembly Tools ---
elif menu == "🧫 Cloning & Assembly Tools":
    run_cloning_tools()

# --- MODULE: Protein Tools ---
elif menu == "🧬 Protein Tools":
    run_protein_tools()

# --- MODULE: Plasmid Designer ---
elif menu == "🧫 Plasmid Designer":
    run_plasmid_designer()

# --- MODULE: Plasmid Plus ---
elif menu == "🧬 Plasmid Plus":
    run_plasmid_plus()

# --- MODULE: Database & Reference Integration ---
elif menu == "🧫 Database & Reference Integration":
    run_database_integration()

# --- MODULE: Data Management ---
elif menu == "🧬 Data Management":
    run_data_management()
    
# --- MODULE: Alignment & Phylogeny ---
elif menu == "🌳 Alignment & Phylogeny":
    run_alignment_phylogeny()


# --- Footer ---
st.markdown("---")
st.caption("🧠 Entwickelt mit ❤️ in Hamburg · Version 2.9 · Zweisprachig DE/EN")
