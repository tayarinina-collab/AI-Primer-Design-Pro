import streamlit as st

# ---- App Configuration ----
st.set_page_config(
    page_title="AI Primer Design Pro",
    page_icon="🧬",
    layout="wide"
)

# ---- Theme Toggle ----
if "theme" not in st.session_state:
    st.session_state["theme"] = "Light"

theme = st.session_state["theme"]

toggle = st.toggle("🌗 Dark / Light Mode", value=(theme == "Dark"))
st.session_state["theme"] = "Dark" if toggle else "Light"

bg_color = "#0D1117" if st.session_state["theme"] == "Dark" else "#F5F7FA"
text_color = "#FFFFFF" if st.session_state["theme"] == "Dark" else "#000000"

st.markdown(
    f"""
    <style>
    body {{
        background-color: {bg_color};
        color: {text_color};
    }}
    .stApp {{
        background-color: {bg_color};
        color: {text_color};
    }}
    </style>
    """,
    unsafe_allow_html=True
)

# ---- Language Selection ----
lang = st.radio("Language / Sprache", ["🇩🇪 Deutsch", "🇬🇧 English"], horizontal=True)

# ---- Titles ----
if lang == "🇩🇪 Deutsch":
    st.title("🧬 Willkommen bei AI Primer Design Pro")
    st.write("""
    Eine intelligente Bioinformatik-Plattform für DNA-, RNA- und Protein-Analysen.  
    Hier vereinen sich KI-gestützte Primer-Entwicklung, Sequenzverwaltung und visuelle Labor-Tools  
    in einer modernen, anpassbaren Umgebung – optimiert für Forschung, Lehre und Innovation.
    """)
else:
    st.title("🧬 Welcome to AI Primer Design Pro")
    st.write("""
    An intelligent bioinformatics platform for DNA, RNA, and protein analysis.  
    Combining AI-driven primer design, sequence management and interactive lab tools  
    in one adaptive, modern environment – built for research, education, and biotech innovation.
    """)

# ---- Grid Overview ----
modules = [
    ("🧬 Sequence Management", "Sequenzverwaltung / Sequence Management"),
    ("🧪 Primer Design & PCR Tools", "Primer-Entwurf & PCR-Simulation"),
    ("⚙️ Batch Processing", "Automatisierte Analyse mehrerer Proben"),
    ("🧫 Cloning & Assembly Tools", "Klonierungs- & Assemblierungs-Assistent"),
    ("💪 Protein Tools", "Protein-Analyse & 3D-Strukturvisualisierung"),
    ("🔗 Database & Reference Integration", "NCBI / UniProt / NEB Verknüpfung"),
    ("🌿 Alignment & Phylogeny", "Sequenzvergleich & phylogenetische Bäume"),
    ("🤖 AI Learning & Chatbot System", "KI-gestützter Labor-Assistent"),
    ("📊 Auto-Report & Visualization", "Automatische Auswertung & Plots"),
    ("🗂️ File Management & Collaboration", "Dateiverwaltung & Team-Freigabe"),
    ("🧠 KI-Innovation & Learning-System", "Adaptives Lern- und Trainingssystem"),
    ("☁️ Cloud Sync & Offline Cache", "Sichere Daten-Synchronisation"),
    ("🔬 Bioinformatics APIs & Integrations", "APIs für BLAST, Primer3, PDB usw."),
    ("⚙️ Settings & User Profiles", "Einstellungen, Themes und Profile")
]

st.markdown("---")

cols = st.columns(4)
for i, (emoji, desc) in enumerate(modules):
    with cols[i % 4]:
        st.markdown(f"### {emoji}")
        st.markdown(desc)
        st.markdown("---")
