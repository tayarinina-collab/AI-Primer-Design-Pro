import streamlit as st

# ------------------ CONFIG ------------------
st.set_page_config(page_title="AI Primer Design Pro", page_icon="🧬", layout="wide")

# ------------------ THEME TOGGLE ------------------
if "theme" not in st.session_state:
    st.session_state["theme"] = "Light"

toggle = st.toggle("🌗 Dark / Light Mode", value=(st.session_state["theme"] == "Dark"))
st.session_state["theme"] = "Dark" if toggle else "Light"

bg_color = "#0D1117" if st.session_state["theme"] == "Dark" else "#F5F7FA"
text_color = "#FFFFFF" if st.session_state["theme"] == "Dark" else "#000000"

st.markdown(
    f"""
    <style>
    body, .stApp {{
        background-color: {bg_color};
        color: {text_color};
        font-family: 'Inter', sans-serif;
    }}
    h1, h2, h3, p {{
        color: {text_color};
    }}
    .center-container {{
        text-align: center;
        padding-top: 40px;
    }}
    .subtitle {{
        font-size: 18px;
        opacity: 0.9;
        line-height: 1.6;
        max-width: 800px;
        margin: auto;
    }}
    </style>
    """,
    unsafe_allow_html=True
)

# ------------------ LANGUAGE TOGGLE ------------------
lang = st.radio("Language / Sprache", ["🇩🇪 Deutsch", "🇬🇧 English"], horizontal=True)

# ------------------ SIDEBAR: MODULE LIST ------------------
with st.sidebar:
    st.markdown("### 🧬 **AI Primer Design Pro**")
    st.write("**Module / Funktionen:**")
    st.markdown("---")

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

    for emoji, desc in modules:
        st.markdown(f"**{emoji}** {desc}")

# ------------------ MAIN CENTER CONTENT ------------------
st.markdown("<div class='center-container'>", unsafe_allow_html=True)

if lang == "🇩🇪 Deutsch":
    st.markdown("## 🧬 Willkommen bei **AI Primer Design Pro**")
    st.markdown(
        """
        <p class='subtitle'>
        Eine intelligente Bioinformatik-Plattform für DNA-, RNA- und Protein-Analysen.<br>
        Hier vereinen sich KI-gestützte Primer-Entwicklung, Sequenzverwaltung und visuelle Labor-Tools<br>
        in einer modernen, anpassbaren Umgebung – optimiert für Forschung, Lehre und Innovation.
        </p>
        """,
        unsafe_allow_html=True
    )
else:
    st.markdown("## 🧬 Welcome to **AI Primer Design Pro**")
    st.markdown(
        """
        <p class='subtitle'>
        An intelligent bioinformatics platform for DNA, RNA, and protein analysis.<br>
        Combining AI-driven primer design, sequence management and interactive lab tools<br>
        in one adaptive, modern environment – built for research, education, and biotech innovation.
        </p>
        """,
        unsafe_allow_html=True
    )

st.markdown("</div>", unsafe_allow_html=True)
