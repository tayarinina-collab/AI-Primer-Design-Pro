# -*- coding: utf-8 -*-
"""
Reports & Export Center
- PDF & CSV Export
- Explainable AI Kommentar
- Zusammenfassung von Modulen
"""

import streamlit as st
import pandas as pd
from io import BytesIO
from reportlab.lib.pagesizes import A4
from reportlab.lib.units import cm
from reportlab.pdfgen import canvas

def run_reports_export_center():
    st.title("📊 Reports & Export Center")
    st.caption("Erstelle vollständige Laborberichte mit KI-Kommentaren und Datenexport")

    # -------------------------------
    # Projektinformationen
    # -------------------------------
    st.markdown("### 🧾 Projektinformationen")
    project_name = st.text_input("Projektname:", "AI Primer Design Pro – Demo Report")
    author = st.text_input("Autor:", "Nina Tayari")
    module = st.selectbox(
        "Datenquelle auswählen:",
        [
            "Primer Design",
            "Protein Tools",
            "Plasmid Designer",
            "Alignment & Phylogeny",
            "Database Integration"
        ]
    )

    st.markdown("---")
    st.markdown("### 🧠 Explainable AI Kommentar")
    ai_comment = st.text_area(
        "Kurze Zusammenfassung oder KI-Begründung:",
        value=f"Das Modul {module} wurde erfolgreich analysiert. Die Daten zeigen stabile Sequenzen, optimale Parameter und konsistente GC-Werte."
    )

    # Beispielhafte DataFrame-Erstellung
    df = pd.DataFrame({
        "Parameter": ["Tm (°C)", "GC %", "Amplicon Länge (bp)", "Bewertung"],
        "Wert": [60, 48, 120, "Sehr gut"]
    })
    st.dataframe(df, use_container_width=True)

    st.markdown("---")
    st.subheader("⬇️ Export Optionen")

    # -------------------------------
    # CSV Export
    # -------------------------------
    csv_bytes = df.to_csv(index=False).encode("utf-8")
    st.download_button(
        label="📄 CSV exportieren",
        data=csv_bytes,
        file_name=f"{project_name.replace(' ', '_')}.csv",
        mime="text/csv"
    )

    # -------------------------------
    # PDF Export
    # -------------------------------
    if st.button("🧾 PDF-Report generieren"):
        pdf_buffer = create_pdf_report(project_name, author, module, ai_comment, df)
        st.download_button(
            label="💾 PDF herunterladen",
            data=pdf_buffer,
            file_name=f"{project_name.replace(' ', '_')}.pdf",
            mime="application/pdf"
        )
        st.success("✅ PDF-Report erfolgreich erstellt!")

    st.markdown("---")
    st.caption("📊 AI Primer Design Pro · Reports & Export Center · Version 1.0")


# -------------------------------
# 📄 PDF-Erstellung
# -------------------------------
def create_pdf_report(project_name, author, module, ai_comment, df):
    buffer = BytesIO()
    c = canvas.Canvas(buffer, pagesize=A4)
    width, height = A4

    # Header
    c.setFont("Helvetica-Bold", 16)
    c.drawString(2*cm, height - 2*cm, "🧬 AI Primer Design Pro – Report")

    # Meta
    c.setFont("Helvetica", 11)
    c.drawString(2*cm, height - 3*cm, f"Projekt: {project_name}")
    c.drawString(2*cm, height - 3.6*cm, f"Autor: {author}")
    c.drawString(2*cm, height - 4.2*cm, f"Modul: {module}")

    # AI Kommentar
    text = c.beginText(2*cm, height - 5.5*cm)
    text.setFont("Helvetica-Oblique", 10)
    text.textLines(f"🧠 KI-Kommentar:\n{ai_comment}")
    c.drawText(text)

    # Tabelle
    y = height - 9*cm
    c.setFont("Helvetica-Bold", 11)
    c.drawString(2*cm, y, "Parameter")
    c.drawString(8*cm, y, "Wert")
    c.setFont("Helvetica", 10)
    y -= 0.4*cm

    for _, row in df.iterrows():
        c.drawString(2*cm, y, str(row["Parameter"]))
        c.drawString(8*cm, y, str(row["Wert"]))
        y -= 0.5*cm

    # Footer
    c.setFont("Helvetica-Oblique", 9)
    c.drawString(2*cm, 2*cm, "Erstellt mit ❤️ in Hamburg · AI Primer Design Pro · Version 2.9")
    c.showPage()
    c.save()
    buffer.seek(0)
    return buffer.getvalue()
