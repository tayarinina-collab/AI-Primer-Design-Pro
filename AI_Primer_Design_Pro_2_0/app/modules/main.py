# ==============================
# 🧬 AI Primer Design Pro – Main App (v2.4 Stable)
# ==============================
import streamlit as st
import importlib
import os, sys

# --- Pfad zu den Modulen hinzufügen ---
sys.path.append(os.path.join(os.path.dirname(__file__), "modules"))

# --- Seitenkonfiguration ---
st.set_page_config(
    page_title="AI Primer Design Pro",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# --- Sprache & Theme Auswahl ---
lang = st.sidebar.radio("🌐 Sprache / Language", ("Deutsch", "English"), horizontal=True)
theme_mode = st.sidebar.radio("🎨 Theme", ("🌙 Dark Mode", "☀️ Light Mode"), horizontal=True)

# --- Dynamische Modul-Erkennung ---
module_dir = os.path.join(os.path.dirname(__file__), "modules")
modules = []

for f in sorted(os.listdir(module_dir)):
    if f.endswith(".py") and f not in ["__init__.py", "main.py"]:
        name = f.replace(".py", "")
        title = name.replace("_", " ").title()
        # optional: Titel aus Kopfzeile # Title:
        try:
            with open(os.path.join(module_dir, f), "r", encoding="utf-8") as mf:
                for line in mf:
                    if line.strip().startswith("# Title:"):
                        title = line.strip().split(":", 1)[1].strip()
                        break
        except:
            pass
        modules.append((title, name))

# --- Sidebar Navigation ---
st.sidebar.markdown("## 🧩 Module auswählen / Select Module")
if not modules:
    st.sidebar.warning("⚠️ Keine Module gefunden. Bitte prüfe den Ordner '/modules'.")
    selected_module = None
else:
    titles = [t[0] for t in modules]
    selected_title = st.sidebar.radio("🔬 Modul wählen / Select Module:", titles)
    selected_module = [m for t, m in modules if t[0] == selected_title][0]

# --- Header ---
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

# --- Theme Styling (funktioniert lokal + in Cloud) ---
if theme_mode == "🌙 Dark Mode":
    st.markdown(
        """
        <style>
        html, body, [data-testid="stAppViewContainer"], [class*="css"], .stApp {
            background-color: #0e1117 !important;
            color: white !important;
        }
        [data-testid="stSidebar"], .stSidebar {
            background-color: #111 !important;
        }
        h1,h2,h3,h4,h5,h6,p,div,span,label {
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
        html, body, [data-testid="stAppViewContainer"], [class*="css"], .stApp {
            background-color: #f8f9fa !important;
            color: #111 !important;
        }
        [data-testid="stSidebar"], .stSidebar {
            background-color: #ffffff !important;
        }
        h1,h2,h3,h4,h5,h6,p,div,span,label {
            color: #111 !important;
        }
        </style>
        """,
        unsafe_allow_html=True,
    )

# --- Modul laden ---
if selected_module:
    try:
        mod = importlib.import_module(f"modules.{selected_module}")
        if hasattr(mod, "render"):
            mod.render()
        else:
            st.warning(f"⚠️ Modul '{selected_module}' gefunden, aber keine render()-Funktion.")
    except Exception as e:
        st.error(f"❌ Fehler beim Laden von '{selected_module}': {e}")

# --- Footer ---
st.markdown(
    """
    <hr>
    <p style='text-align:center; color:gray; font-size:14px;'>
    🧠 Developed with ❤️ in Hamburg · Version 2.4 · Bilingual DE/EN
    </p>
    """,
    unsafe_allow_html=True,
)
