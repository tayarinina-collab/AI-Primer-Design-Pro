import streamlit as st
import random

def set_theme():
    # Eindeutiger Key, um Streamlit-Fehler zu vermeiden
    unique_key = f"theme_toggle_{random.randint(1000,9999)}"

    # Sidebar Theme Toggle
    st.sidebar.markdown("## 🧪 Theme & Sprache")
    dark_mode = st.sidebar.toggle("🌗 Dark / Light Mode", key=unique_key)

    # Theme speichern
    if dark_mode:
        st.session_state["theme"] = "dark"
        st.markdown(
            """
            <style>
            body {background-color: #0E1117; color: white;}
            .stButton button {background-color: #1E1E1E; color: white;}
            </style>
            """, unsafe_allow_html=True)
    else:
        st.session_state["theme"] = "light"
        st.markdown(
            """
            <style>
            body {background-color: white; color: black;}
            .stButton button {background-color: #E0E0E0; color: black;}
            </style>
            """, unsafe_allow_html=True)

def set_language():
    # Sprachwahl (Deutsch/Englisch)
    st.sidebar.markdown("## 🌍 Sprache / Language")
    lang = st.sidebar.radio("Wähle Sprache:", ["Deutsch", "English"], key=f"lang_{random.randint(1000,9999)}")
    st.session_state["lang"] = "DE" if lang == "Deutsch" else "EN"
