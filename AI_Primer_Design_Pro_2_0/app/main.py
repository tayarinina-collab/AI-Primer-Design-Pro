import streamlit as st

# ---------- App Config ----------
st.set_page_config(page_title="AI Primer Design Pro", page_icon="🧬", layout="wide")

# ---------- Sidebar ----------
with st.sidebar:
    st.markdown("### 🧬 AI Primer Design Pro")

    # Theme toggle inside sidebar
    if "theme" not in st.session_state:
        st.session_state["theme"] = "Light"

    toggle = st.toggle("🌗 Dark / Light Mode", value=(st.session_state["theme"] == "Dark"))
    st.session_state["theme"] = "Dark" if toggle else "Light"

    # Sidebar module list (numbered)
    st.markdown("---")
    st.markdown("### ⚙️ Funktionen / Modules")

    modules = [
        "1️⃣ Sequence Management",
        "2️⃣ Primer Design & PCR Tools",
        "3️⃣ Batch Processing",
        "4️⃣ Cloning & Assembly Tools",
        "5️⃣ Protein Tools",
        "6️⃣ Database & Reference Integration",
        "7️⃣ Alignment & Phylogeny",
        "8️⃣ AI Learning & Chatbot System",
        "9️⃣ Auto-Report & Visualization",
        "🔟 File Management & Collaboration",
        "11️⃣ KI-Innovation & Learning-System",
        "12️⃣ Cloud Sync & Offline Cache",
        "13️⃣ Bioinformatics APIs & Integrations",
        "14️⃣ Settings & User Profiles"
    ]

    for item in modules:
        st.markdown(f"- {item}")

# ---------- Theme Colors ----------
if st.session_state["theme"] == "Dark":
    bg_color = "#0D1117"
    text_color = "#FFFFFF"
else:
    bg_color = "#F5F7FA"
    text_color = "#000000"

# ---------- Custom CSS ----------
st.markdown(
    f"""
    <style>
    body, .stApp {{
        background-color: {bg_color};
        color: {text_color};
        font-family: 'Inter', sans-serif;
    }}
    h1, h2, h3, p, li, label {{
        color: {text_color};
    }}
    .center-container {{
        text-align: center;
        padding-top: 60px;
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

# ---------- Language Toggle ----------
lang = st.radio("Language / Sprache", ["🇩🇪 Deutsch", "🇬🇧 English"], horizontal=True)

# ---------- Main Center Content ----------
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
