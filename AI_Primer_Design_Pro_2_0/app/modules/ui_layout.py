import streamlit as st

def set_theme():
    # Sidebar Theme Toggle (mit festem Key)
    st.sidebar.markdown("## 🧪 Theme & Sprache")
    dark_mode = st.sidebar.toggle("🌗 Dark / Light Mode", key="theme_toggle")

    # Speichere Zustand in Session State
    if "theme" not in st.session_state:
        st.session_state["theme"] = "light"

    st.session_state["theme"] = "dark" if dark_mode else "light"

    # Wende Theme-Stile an
    if st.session_state["theme"] == "dark":
        st.markdown(
            """
            <style>
            .stApp, .stMarkdown, div, p, h1, h2, h3, h4, h5, h6, label, span {
                color: white !important;
            }
            .stApp {
                background-color: #0E1117 !important;
            }
            .stSidebar {
                background-color: #111111 !important;
            }
            .stButton>button {
                background-color: #1E1E1E !important;
                color: white !important;
                border: 1px solid #333 !important;
            }
            </style>
            """,
            unsafe_allow_html=True
        )
    else:
        st.markdown(
            """
            <style>
            .stApp, .stMarkdown, div, p, h1, h2, h3, h4, h5, h6, label, span {
                color: #111111 !important;
            }
            .stApp {
                background-color: #F8F9FA !important;
            }
            .stSidebar {
                background-color: #FFFFFF !important;
            }
            .stButton>button {
                background-color: #E0E0E0 !important;
                color: #111 !important;
                border: 1px solid #CCC !important;
            }
            </style>
            """,
            unsafe_allow_html=True
        )

def set_language():
    # Sprachwahl (Deutsch/Englisch)
    st.sidebar.markdown("## 🌍 Sprache / Language")
    lang = st.sidebar.radio("Wähle Sprache:", ["Deutsch", "English"], key="lang_toggle")
    st.session_state["lang"] = "DE" if lang == "Deutsch" else "EN"
