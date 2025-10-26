# ==============================
# 🧬 AI Primer Design Pro – Main App (v2.2)
# ==============================
import streamlit as st
import importlib, os, sys

# --- Pfad hinzufügen ---
sys.path.append(os.path.join(os.path.dirname(__file__), "modules"))

# --- Page config ---
st.set_page_config(page_title="AI Primer Design Pro", page_icon="🧬", layout="wide")

# --- Sidebar: Sprache & Theme ---
lang = st.sidebar.radio("🌐 Sprache / Language", ("Deutsch", "English"), horizontal=True)
theme_mode = st.sidebar.radio("🎨 Theme", ("🌙 Dark Mode", "☀️ Light Mode"), horizontal=True)

# --- Theme-Styles dynamisch ---
if theme_mode == "🌙 Dark Mode":
    bg, text = "#0e1117", "white"
else:
    bg, text = "#f8f9fa", "#111"

st.markdown(
    f"""
    <style>
    .stApp, body {{
        background-color: {bg} !important;
        color: {text} !important;
    }}
    h1, h2, h3, h4, h5, h6, p, div, span, label {{
        color: {text} !important;
    }}
    </style>
    """,
    unsafe_allow_html=True,
)

# --- Dynamische Modul-Erkennung ---
module_dir = os.path.join(os.path.dirname(__file__), "modules")
available_modules = []

for f in sorted(os.listdir(module_dir)):
    if f.endswith(".py") and f not in ["__init__.py"]:
        file_path = os.path.join(module_dir, f)
        name = f[:-3]
        # Titel automatisch aus Header-Zeile ziehen
        title = f.replace("_", " ").replace(".py", "").title()
        with open(file_path, "r", encoding="utf-8") as mf:
            for line in mf:
                if line.strip().startswith("# Title:"):
                    title = line.strip().split(":", 1)[1].strip()
                    break
        available_modules.append((title, name))

# --- Sidebar Auswahl ---
st.sidebar.markdown("## 🧩 Module")
titles = [t[0] for t in available_modules]
selected_title = st.sidebar.radio("🔬 Modul wählen:", titles)
selected_module = [m for t, m in available_modules if t == selected_title][0]

# --- Header ---
st.markdown(
    f"""
    <h1 style='text-align:center;'>🧬 AI Primer Design Pro</h1>
    <p style='text-align:center; font-size:18px;'>
        Intelligente Bioinformatik-Plattform für DNA-, RNA- & Protein-Analysen.<br>
        Combining AI, Automation & Visualization for smarter research.
    </p>
    """,
    unsafe_allow_html=True,
)

# --- Modul laden ---
try:
    mod = importlib.import_module(f"modules.{selected_module}")
    if hasattr(mod, "render"):
        mod.render()
    else:
        st.warning(f"⚠️ Modul '{selected_module}' gefunden, aber keine render()-Funktion vorhanden.")
except Exception as e:
    st.error(f"❌ Fehler beim Laden von '{selected_module}': {e}")

# --- Footer ---
st.markdown(
    """
    <hr>
    <p style='text-align:center; color:gray; font-size:14px;'>
    🧠 Developed with ❤️ in Hamburg · Version 2.2 · Bilingual DE/EN
    </p>
    """,
    unsafe_allow_html=True,
)
