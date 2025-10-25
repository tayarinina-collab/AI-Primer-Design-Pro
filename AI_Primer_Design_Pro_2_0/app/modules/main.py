# ==============================
# AI Primer Design Pro – Main App
# ==============================
import streamlit as st
import importlib
import sys, os

# --- Set working path so modules can be found ---
sys.path.append(os.path.join(os.path.dirname(__file__), "modules"))

# --- Page configuration ---
st.set_page_config(
    page_title="AI Primer Design Pro",
    layout="wide",
    page_icon="🧬",
)

# --- Language selection ---
lang = st.sidebar.radio("🌐 Language / Sprache", ("Deutsch", "English"))

# --- Theme switch ---
theme_mode = st.sidebar.radio("🎨 Theme", ("Dark Mode", "Light Mode"))
if theme_mode == "Dark Mode":
    st.markdown(
        """
        <style>
        body { background-color: #0e1117; color: white; }
        .stApp { background-color: #0e1117; color: white; }
        </style>
        """,
        unsafe_allow_html=True,
    )

# --- App title ---
st.markdown(
    """
    <h1 style='text-align:center;'>🧬 AI Primer Design Pro</h1>
    <p style='text-align:center; font-size:18px;'>
    Intelligente Bioinformatik-Plattform für Sequenzanalyse, Primer-Design & Visualisierung.
    </p>
    """,
    unsafe_allow_html=True,
)

# --- Sidebar Navigation ---
st.sidebar.markdown("## 🔬 Module")

modules = {
    "Sequence Management": "sequence_management",
    "Primer Design": "primer_design",
    "Protein Tools": "protein_tools",
    "Phylogeny": "phylogeny",
    "Alignments": "alignments",
    "Restriction Tools": "restriction_tools",
    "Reports": "reports",
    "Settings / About": "settings_about",
}

choice = st.sidebar.radio("🧩 Select Module", list(modules.keys()))

# --- Load selected module dynamically ---
try:
    selected_module = modules[choice]
    module = importlib.import_module(f"modules.{selected_module}")

    # Each module must contain a function run_<modulename>()
    run_function_name = f"run_{selected_module}"
    if hasattr(module, run_function_name):
        getattr(module, run_function_name)()
    else:
        st.warning(f"⚠️ Modul '{choice}' wurde gefunden, aber keine Funktion '{run_function_name}()' existiert.")
except Exception as e:
    st.error(f"❌ Fehler beim Laden des Moduls '{choice}': {e}")

# --- Footer ---
st.markdown(
    """
    <hr>
    <p style='text-align:center; color:gray; font-size:14px;'>
    v2.0 Beta · Zweisprachig DE/EN · Entwickelt mit Streamlit
    </p>
    """,
    unsafe_allow_html=True,
)
