import streamlit as st
import pandas as pd
from io import BytesIO
from reportlab.lib.pagesizes import A4
from reportlab.lib.units import cm
from reportlab.pdfgen import canvas

def run_reports_export_center():
    st.title("📊 Reports & Export Center")
    st.caption("Erstelle vollständige Laborberichte mit KI-Kommentaren und Datenexport")

    st.markdown("### 🧾 Projektinformationen")
    project_name = st.text_input("Projektname:", "AI Primer Design Pro – Demo Report")
    author = st.text_input("Autor:", "Nina Tayari")
    module = st.selectbox(
        "Datenquelle auswählen:",
        ["Primer Design", "Protein Tools", "Plasmid Designer", "Alignment & Phylogeny", "Database Integration"]
    )

    st.markdown("---")
    ai_comment = st.text_area(
        "🧠 Explainable AI Kommentar:",
        value=f"Das Modul {module} wurde erfolgreich analysiert. Die Daten zeigen stabile Sequenzen und optimale Parameter."
    )

    df = pd.DataFrame({
        "Parameter": ["Tm (°C)", "GC %", "Amplicon Länge (bp)", "Bewertung"],
        "Wert": [60, 48, 120, "Sehr gut"]
    })
    st.dataframe(df, use_container_width=True)

    csv_bytes = df.to_csv(index=False).encode("utf-8")
    st.download_button("📄 CSV exportieren", csv_bytes, f"{project_name.replace(' ', '_')}.csv", "text/csv")

    if st.button("🧾 PDF-Report generieren"):
        pdf_data = create_pdf_report(project_name, author, module, ai_comment, df)
        st.download_button("💾 PDF herunterladen", pdf_data, f"{project_name.replace(' ', '_')}.pdf", "application/pdf")
        st.success("✅ PDF-Report erfolgreich erstellt!")

def create_pdf_report(project_name, author, module, ai_comment, df):
    buffer = BytesIO()
    c = canvas.Canvas(buffer, pagesize=A4)
    width, height = A4

    c.setFont("Helvetica-Bold", 16)
    c.drawString(2*cm, height - 2*cm, "🧬 AI Primer Design Pro – Report")
    c.setFont("Helvetica", 11)
    c.drawString(2*cm, height - 3*cm, f"Projekt: {project_name}")
    c.drawString(2*cm, height - 3.6*cm, f"Autor: {author}")
    c.drawString(2*cm, height - 4.2*cm, f"Modul: {module}")

    text = c.beginText(2*cm, height - 5.5*cm)
    text.setFont("Helvetica-Oblique", 10)
    text.textLines(f"🧠 KI-Kommentar:\n{ai_comment}")
    c.drawText(text)

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

    c.setFont("Helvetica-Oblique", 9)
    c.drawString(2*cm, 2*cm, "Erstellt mit ❤️ in Hamburg · AI Primer Design Pro · Version 2.9")
    c.showPage()
    c.save()
    buffer.seek(0)
    return buffer.getvalue()
