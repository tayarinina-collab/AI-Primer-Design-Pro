import streamlit as st

# ---------- CONFIG ----------
st.set_page_config(page_title="AI Primer Design Pro", page_icon="🧬", layout="wide")

# ---------- SIDEBAR ----------
with st.sidebar:
    st.markdown("### 🧬 **AI Primer Design Pro**")

    # Theme toggle inside sidebar
    if "theme" not in st.session_state:
        st.session_state["theme"] = "Light"

    toggle = st.toggle("🌗 Dark / Light Mode", value=(st.session_state["theme"] == "Dark"))
    st.session_state["theme"] = "Dark" if toggle else "Light"

    # Module list (with emojis)
    st.markdown("---")
    st.markdown("### 🧪 Funktionen / Modules")

    modules = [
        "🧬 Sequence Management",
        "🧪 Primer Design & PCR Tools",
        "⚙️ Batch Processing",
        "🧫 Cloning & Assembly Tools",
        "💪 Protein Tools",
        "🔗 Database & Reference Integration",
        "🌿 Alignment & Phylogeny",
        "🤖 AI Learning & Chatbot System",
        "📊 Auto-Report & Visualization",
        "🗂️ File Management & Collaboration",
        "🧠 KI-Innovation & Learning-System",
        "☁️ Cloud Sync & Offline Cache",
        "🔬 Bioinformatics APIs & Integrations",
        "⚙️ Settings & User Profiles"
    ]

    for item in modules:
        st.markdown(f"- {item}")

# ---------- COLORS ----------
if st.session_state["theme"] == "Dark":
    bg_color = "#0D1117"        # Dark background
    sidebar_color = "#161B22"   # Dark sidebar
    text_color = "#FFFFFF"
else:
    bg_color = "#F5F7FA"        # Light background
    sidebar_color = "#FFFFFF"   # Light sidebar
    text_color = "#000000"

# ---------- CUSTOM CSS ----------
st.markdown(
    f"""
    <style>
    /* General background + text */
    body, .stApp {{
        background-color: {bg_color};
        color: {text_color};
        font-family: 'Inter', sans-serif;
    }}
    h1, h2, h3, p, li, label {{
        color: {text_color};
    }}

    /* Sidebar styling */
    section[data-testid="stSidebar"] {{
        background-color: {sidebar_color};
        color: {text_color};
        border-right: 1px solid rgba(255,255,255,0.1);
        padding-left: 10px;
    }}
    .stSidebar .stMarkdown, .stSidebar p, .stSidebar li, .stSidebar label {{
        color: {text_color} !important;
    }}

    /* Center alignment for main content */
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

# ---------- LANGUAGE TOGGLE ----------
lang = st.radio("Language / Sprache", ["🇩🇪 Deutsch", "🇬🇧 English"], horizontal=True)

# ---------- CENTERED INTRO TEXT ----------
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
