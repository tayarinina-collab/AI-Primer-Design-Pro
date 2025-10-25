import streamlit as st
from pathlib import Path

st.set_page_config(page_title="AI Primer Design Pro", page_icon="🧬", layout="wide")
styles_path = Path(__file__).parent / "assets" / "styles.css"
if styles_path.exists():
    st.markdown("<style>" + styles_path.read_text() + "</style>", unsafe_allow_html=True)

if "lang" not in st.session_state:
    st.session_state["lang"] = "de"  # default: German for presentation

def t(key: str) -> str:
    L = {
        "title": {"en": "AI Primer Design Pro", "de": "AI Primer Design Pro"},
        "subtitle": {
            "en": "Intelligent bioinformatics platform for sequence analysis, primer design & visualization.",
            "de": "Intelligente Bioinformatik-Plattform fuer Sequenzanalyse, Primer-Design & Visualisierung.",
        },
        "nav": {
            "en": ["Primer Design", "In-Silico PCR", "Plasmid Map", "Restriction Tools",
                   "Protein Tools", "Alignments", "Phylogeny", "Reports", "Settings / About"],
            "de": ["Primer-Design", "In-Silico PCR", "Plasmid-Karte", "Restriktions-Tools",
                   "Protein-Tools", "Alignments", "Phylogenie", "Reports", "Einstellungen / Info"],
        },
        "footer": {
            "en": "v2.0 Beta · Built with Streamlit · Bilingual EN/DE",
            "de": "v2.0 Beta · Entwickelt mit Streamlit · Zweisprachig DE/EN",
        },
    }
    return L.get(key, {}).get(st.session_state["lang"], key)

with st.sidebar:
    st.markdown("### 🌐 Language / Sprache")
    lang = st.radio("Language", ["de", "en"], format_func=lambda x: "Deutsch 🇩🇪" if x=="de" else "English 🇬🇧", index=0)
    st.session_state["lang"] = lang
    st.markdown("---")
    st.markdown("**Modules / Module**")
    st.caption("Select a tab from the top navigation.")
    st.markdown("---")
    st.caption(t("footer"))

st.title("🧬 " + t("title"))
st.write(t("subtitle"))

tabs = st.tabs(t("nav"))

from modules import primer_design, in_silico_pcr, plasmid_designer, restriction_tools
from modules import protein_tools, alignments, phylogeny, reports, settings_about

with tabs[0]:
    primer_design.render()

with tabs[1]:
    in_silico_pcr.render()

with tabs[2]:
    plasmid_designer.render()

with tabs[3]:
    restriction_tools.render()

with tabs[4]:
    protein_tools.render()

with tabs[5]:
    alignments.render()

with tabs[6]:
    phylogeny.render()

with tabs[7]:
    reports.render()

with tabs[8]:
    settings_about.render()