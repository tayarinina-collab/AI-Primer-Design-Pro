# ==============================
# 🧬 AI Primer Design Pro – Main App (v2.1)
# ==============================
import streamlit as st
import importlib
import sys, os

# --- Make sure the 'modules' path is available ---
sys.path.append(os.path.join(os.path.dirname(__file__), "modules"))

# --- Page config ---
st.set_page_config(
    page_title="AI Primer Design Pro",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# --- Sidebar: Sprache / Language ---
lang = st.sidebar.radio("🌐 Sprache / Language", ("Deutsch", "English"), horizontal=True)

# --- Sidebar: Theme Switch ---
theme_mode = st.sidebar.radio("🎨 Theme", ("🌙 Dark Mode", "☀️ Light Mode"), horizontal=True)

# --- Sidebar Navigation ---
st.sidebar.markdown("## 🧩 Module")

modules = {
    "🏠 Overview / Übersicht": "overview",
    "🧬 Sequence Management": "sequence_management",
    "🧫 Primer Design": "primer_design",
    "🧪 Primer Design – Advanced": "primer_design_advanced",
    "🧫 Cloning & Assembly Tools": "cloning_tools",
    "🧬 Protein Tools": "protein_tools",
    "🧫 Database & Reference Integration": "database_integration",
    "🧫 Plasmid Karte": "plasmid_designer",
    "🧬 Plasmid Plus": "plasmid_plus",
    "📊 Reports": "reports",
    "⚙️ Settings / About": "settings_about",
}

choice = st.sidebar.radio("🔬 Select Module", list(modules.keys()))

# --- Theme CSS ---
if theme_mode == "🌙 Dark Mode":
    st.markdown(
        """
        <style>
        .stApp, body {
            background-color: #0e1117 !important;
            color: white !important;
        }
        .stSidebar {
            background-color: #111 !important;
        }
        h1, h2, h3, h4, h5, h6, p, div, span, label {
            color: white !important;
        }
        </style>
        """,
        unsafe_allow_html=True,
    )
else:
    st.markdown(
        """
        <style>
        .stApp, body {
            background-color: #f8f9fa !important;
            color: #111 !important;
        }
        .stSidebar {
            background-color: #ffffff !important;
        }
        h1, h2, h3, h4, h5, h6, p, div, span, label {
            color: #111 !important;
        }
        </style>
        """,
        unsafe_allow_html=True,
    )

# --- App Header ---
st.markdown(
    """
    <h1 style='text-align:center;'>🧬 AI Primer Design Pro</h1>
    <p style='text-align:center; font-size:18px;'>
        Intelligente Bioinformatik-Plattform für DNA-, RNA- & Protein-Analysen.<br>
        Combining AI, Automation & Visualization for smarter research.
    </p>
    """,
    unsafe_allow_html=True,
)

# --- Dynamische Modul-Ladung ---
try:
    selected_module = modules[choice]
    module = importlib.import_module(f"modules.{selected_module}")

    run_function_name = f"run_{selected_module}"
    if hasattr(module, run_function_name):
        getattr(module, run_function_name)()
    else:
        st.warning(f"⚠️ Modul '{choice}' gefunden, aber keine Funktion '{run_function_name}()' in der Datei.")
except Exception as e:
    st.error(f"❌ Fehler beim Laden von '{choice}': {e}")

# --- Footer ---
st.markdown(
    """
    <hr>
    <p style='text-align:center; color:gray; font-size:14px;'>
    🧠 Developed with ❤️ in Hamburg · Version 2.1 · Bilingual DE/EN
    </p>
    """,
    unsafe_allow_html=True,
)
