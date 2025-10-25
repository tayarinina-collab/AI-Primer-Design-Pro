import streamlit as st

def set_theme():
    if "theme" not in st.session_state:
        st.session_state["theme"] = "Light"

    toggle = st.sidebar.toggle("🌗 Dark / Light Mode", value=(st.session_state["theme"] == "Dark"))
    st.session_state["theme"] = "Dark" if toggle else "Light"

    bg_color = "#0D1117" if st.session_state["theme"] == "Dark" else "#F5F7FA"
    text_color = "#FFFFFF" if st.session_state["theme"] == "Dark" else "#000000"

    st.markdown(
        f"""
        <style>
        body, .stApp {{ background-color: {bg_color}; color: {text_color}; }}
        h1, h2, h3, p, li, label {{ color: {text_color}; }}
        section[data-testid="stSidebar"] {{ background-color: {bg_color}; }}
        </style>
        """,
        unsafe_allow_html=True
    )
