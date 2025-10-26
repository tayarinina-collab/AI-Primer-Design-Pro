# ==============================
# 🧬 AI Primer Design Pro – Main App (v2.5 Debug)
# ==============================
import streamlit as st
import importlib
import os, sys, traceback

# --- Pfad zu Modulen hinzufügen ---
sys.path.append(os.path.join(os.path.dirname(__file__), "modules"))

# --- Seitenkonfiguration ---
st.set_page_config(
    page_title="AI Primer Design Pro",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# --- Sprache & Theme ---
lang = st.sidebar.radio("🌐 Sprache / Language", ("Deutsch", "English"), horizontal=True)
theme_mode = st.sidebar.radio("🎨 Theme", ("🌙 Dark Mode", "☀️ Light Mode"), horizontal=True)

# --- Theme CSS ---
if theme_mode == "🌙 Dark Mode":
    st.markdown("""
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
    """, unsafe_allow_html=True)
else:
    st.markdown("""
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
    """, unsafe_allow_html=True)

# --- Module automatisch erkennen ---
module_dir = os.path.join(os.path.dirname(__file__), "modules")
modules_found = []
modules_errors = {}

for f in sorted(os.listdir(module_dir)):
    if f.endswith(".py") and f not in ["__init__.py", "main.py"]:
        name = f.replace(".py", "")
        title = name.replace("_", " ").title()
        try:
            # Optional: Titel aus Kopfkommentar
            with open(os.path.join(module_dir, f), "r", encoding="utf-8") as mf:
                for line in mf:
                    if line.strip().startswith("# Title:"):
                        title = line.strip().split(":", 1)[1].strip()
                        break

            # Versuch zu importieren
            importlib.import_module(f"modules.{name}")
            modules_found.append((title, name))
        except Exception as e:
            modules_errors[name] = str(e)
            modules_found.append((f"❌ {title} (Fehler beim Import)", None))

# --- Sidebar Navigation ---
st.sidebar.markdown("## 🧩 Module auswählen / Select Module")
valid_titles = [t for t, n in modules_found if n is not None]
invalid_titles = [t for t, n in modules_found if n is None]

if not valid_titles and not invalid_titles:
    st.sidebar.warning("⚠️ Keine Module gefunden – prüfe den Ordner 'modules/'.")
else:
    selected_title = st.sidebar.radio("🔬 Modul wählen / Select Module:", [t for t, _ in modules_found])
    selected_entry = [m for t, m in modules_found if t == selected_title]
    selected_module = selected_entry[0] if selected_entry else None

# --- Header ---
st.markdown("""
<h1 style='text-align:center;'>🧬 AI Primer Design Pro</h1>
<p style='text-align:center; font-size:18px;'>
    Intelligente Bioinformatik-Plattform für DNA-, RNA- & Protein-Analysen.<br>
    Combining AI, Automation & Visualization for smarter research.
</p>
""", unsafe_allow_html=True)

# --- Modul laden ---
if selected_module:
    try:
        mod = importlib.import_module(f"modules.{selected_module}")
        if hasattr(mod, "render"):
            mod.render()
        else:
            st.warning(f"⚠️ Modul '{selected_module}' gefunden, aber keine render()-Funktion.")
    except Exception as e:
        st.error(f"❌ Fehler beim Ausführen von '{selected_module}': {e}")
        st.exception(e)
else:
    if "Fehler" in selected_title:
        st.error("Dieses Modul konnte nicht geladen werden – Details unten im Debug-Bereich.")

# --- Debug-Bereich unten ---
st.markdown("---")
st.subheader("🧩 Modul-Diagnose")
if modules_errors:
    for mod_name, err in modules_errors.items():
        st.error(f"❌ **{mod_name}.py** konnte nicht geladen werden:\n```\n{err}\n```")
else:
    st.success("✅ Alle Module wurden erfolgreich importiert!")

# --- Footer ---
st.markdown("""
<hr>
<p style='text-align:center; color:gray; font-size:14px;'>
🧠 Developed with ❤️ in Hamburg · Version 2.5 · Debug Mode Enabled
</p>
""", unsafe_allow_html=True)
