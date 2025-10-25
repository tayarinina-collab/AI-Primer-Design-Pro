import streamlit as st
from modules.sequence_tools import compute_basic_properties, find_orfs, find_motifs, gc_profile
from modules.ai_assistant import interpret_sequence
from modules.ui_layout import set_theme
import matplotlib.pyplot as plt

def run_sequence_management():
    set_theme()
    lang = st.radio("Language / Sprache", ["🇩🇪 Deutsch", "🇬🇧 English"], horizontal=True)

    st.title("🧬 Sequence Management")
    st.markdown("Upload or paste your sequence below:" if lang == "🇬🇧 English" else "Lade deine Sequenz hoch oder füge sie unten ein:")

    uploaded_file = st.file_uploader("Upload File", type=["fasta", "gb", "txt"])
    seq_input = st.text_area("🧫 Sequence Input")

    if uploaded_file:
        seq = uploaded_file.read().decode("utf-8")
    else:
        seq = seq_input.strip()

    if seq:
        st.markdown("---")
        st.subheader("🧩 Sequence Analysis / Sequenzanalyse")

        props = compute_basic_properties(seq)
        st.write(props)

        orfs = find_orfs(seq)
        motifs = find_motifs(seq)
        st.write("📍 ORFs:", orfs)
        st.write("🔎 Motifs:", motifs)

        # GC Profile Plot
        df = gc_profile(seq)
        fig, ax = plt.subplots()
        ax.plot(df["Window"], df["GC%"])
        ax.set_xlabel("Window")
        ax.set_ylabel("GC%")
        ax.set_title("GC Profile")
        st.pyplot(fig)

        # KI-Assistent
        st.markdown("### 🤖 KI-Assistent / AI Assistant")
        with st.spinner("Analysiere Sequenz..."):
            result = interpret_sequence(seq, lang="DE" if lang == "🇩🇪 Deutsch" else "EN")
        st.success(result)

        # Export
        st.download_button("⬇️ Export as FASTA", data=seq, file_name="sequence.fasta")
