# -*- coding: utf-8 -*-
import streamlit as st
import pandas as pd
from io import StringIO

def run_data_management():
    st.title("🧬 Data Management")
    st.caption("Verwaltung, Upload und Export biologischer Daten")

    st.markdown("### 📤 Datei-Upload")
    uploaded_file = st.file_uploader(
        "Datei hochladen (CSV, FASTA oder GenBank)",
        type=["csv", "txt", "fasta", "fa", "gb"]
    )

    if uploaded_file:
        file_type = uploaded_file.name.split(".")[-1].lower()

        if file_type == "csv":
            df = pd.read_csv(uploaded_file)
            st.success("✅ CSV-Datei erfolgreich geladen!")
            st.dataframe(df, use_container_width=True)
            csv_export = df.to_csv(index=False).encode("utf-8")
            st.download_button(
                "💾 CSV exportieren",
                data=csv_export,
                file_name="data_export.csv",
                mime="text/csv"
            )

        elif file_type in ["fasta", "fa", "txt"]:
            text = uploaded_file.read().decode("utf-8")
            st.success("✅ FASTA-Datei erfolgreich geladen!")
            st.text_area("📄 Inhalt", text, height=200)

        elif file_type == "gb":
            text = uploaded_file.read().decode("utf-8")
            st.success("✅ GenBank-Datei erfolgreich geladen!")
            st.text_area("📄 Inhalt", text, height=200)

        else:
            st.warning("⚠️ Unbekanntes Dateiformat. Bitte CSV, FASTA oder GenBank hochladen.")
    else:
        st.info("Bitte eine Datei hochladen, um sie hier zu verwalten.")

    st.markdown("---")
    st.markdown("### 📊 Vorschau / Export")
    st.caption("Hier können Daten künftig zwischen Modulen geteilt und synchronisiert werden.")
