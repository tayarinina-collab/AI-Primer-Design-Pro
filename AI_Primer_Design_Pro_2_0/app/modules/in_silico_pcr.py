import streamlit as st

def run_in_silico_pcr():
    st.title("🧫 In-Silico PCR Modul")
    st.write("Hier kannst du virtuelle PCR-Simulationen durchführen.")

    uploaded_file = st.file_uploader("🔬 Lade eine DNA-Sequenzdatei hoch (FASTA/GenBank):", type=["fasta", "gb", "txt"])
    
    if uploaded_file:
        content = uploaded_file.read().decode("utf-8")
        st.text_area("📄 Eingelesene Sequenz:", content, height=200)
        st.success("✅ Datei erfolgreich geladen.")
    else:
        st.info("Bitte lade eine Sequenzdatei hoch, um fortzufahren.")
