# -*- coding: utf-8 -*-
import streamlit as st

# --- Module Imports ---
from modules.reports_export_center import run_reports_export_center
from modules.ai_learning_chatbot import run_ai_learning_chatbot
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
from modules.gene_map_viewer import visualize_dna_map   # 🧬 NEW IMPORT

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
        "🧬 Visual DNA Map",                  # 🧩 NEW MODULE
        "🧫 Cloning & Assembly Tools",
        "🧬 Protein Tools",
        "🧫 Plasmid Designer",
        "🧬 Plasmid Plus",
        "🧫 Database & Reference Integration",
        "🧬 Data Management",
        "🌳 Alignment & Phylogeny",
        "🤖 AI Learning & Chatbot System",
        "📊 Reports & Export Center",
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

# --- MODULE: Visual DNA Map ---
elif menu == "🧬 Visual DNA Map":
    st.title("🧬 Visual DNA Map & Primer Heatmap")

    fasta_file = st.file_uploader("📂 Lade eine FASTA-Datei hoch", type=["fasta", "fa"])
    if fasta_file:
        with open("uploaded.fasta", "wb") as f:
            f.write(fasta_file.getbuffer())

        # Beispiel-Primer (später dynamisch aus Primer-Design übernehmen)
        primers = [
            {"name": "Fwd1", "start": 120, "end": 140, "Tm": 59.2, "GC": 45, "score": 90},
            {"name": "Rev1", "start": 420, "end": 440, "Tm": 61.5, "GC": 52, "score": 70},
        ]

        st.success("✅ Datei geladen! DNA-Karte wird generiert...")
        visualize_dna_map("uploaded.fasta", primers, color_by="score")

    else:
        st.info("⬆️ Bitte lade eine FASTA-Datei hoch, um die DNA-Karte anzuzeigen.")

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

# --- MODULE: AI Learning & Chatbot System ---
elif menu == "🤖 AI Learning & Chatbot System":
    run_ai_learning_chatbot()

# --- MODULE: Reports & Export Center ---
elif menu == "📊 Reports & Export Center":
    run_reports_export_center()

# --- Footer ---
st.markdown("---")
st.caption("🧠 Entwickelt mit ❤️ in Hamburg · Version 2.9 · Zweisprachig DE/EN")
