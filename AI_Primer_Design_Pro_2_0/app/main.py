# -*- coding: utf-8 -*-
import streamlit as st

# --- Module Imports ---
from modules.reports_export_center import run_reports_export_center
from modules.ai_learning_chatbot import run_ai_learning_chatbot
from modules.alignment_phylogeny import run_alignment_phylogeny
from modules.sequence_management import run_sequence_management
from modules.primer_design import run_primer_design
from modules.primer_design_advanced import run_primer_design_advanced
from modules.cloning_tools import run_cloning_tools
from modules.protein_tools import run_protein_tools
from modules.plasmid_designer import run_plasmid_designer
from modules.plasmid_plus import run_plasmid_plus
from modules.database_integration import run_database_integration
from modules.data_management import run_data_management
from modules.ui_layout import set_theme

# --- Seiteneinstellungen ---
st.set_page_config(
    page_title="AI Primer Design Pro",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# --- Theme Setup ---
set_theme()

# --- Sidebar Titel ---
st.sidebar.markdown("<h2 style='text-align:center;'>🧬 AI Primer Design Pro</h2>", unsafe_allow_html=True)
st.sidebar.caption("Intelligente Bioinformatik-Plattform für moderne Labore")
st.sidebar.markdown("---")

# --- Sprachumschalter ---
language = st.sidebar.radio("🌍 Sprache / Language", ["Deutsch", "English"], horizontal=True)
st.sidebar.markdown("---")

# --- Navigation mit visuellen Icons ---
st.sidebar.markdown("### 🧩 Module auswählen / Select Module")
menu = st.sidebar.radio(
    "Modul-Liste",
    [
        "🏠 Übersicht",
        "🧬 Sequence Management",
        "🧫 Primer Design",
        "🧪 Primer Design – Advanced",
        "🧬 Visual DNA Map",
        "🧫 Cloning & Assembly Tools",
        "🧬 Protein Tools",
        "🧫 Plasmid Designer",
        "🧬 Plasmid Plus",
        "🧫 Database & Reference Integration",
        "🧬 Data Management",
        "🌳 Alignment & Phylogeny",
        "🤖 AI Learning & Chatbot System",
        "📊 Reports & Export Center",
    ],
)

# --- HAUPTINHALT ---
if menu == "🏠 Übersicht":
    if language == "Deutsch":
        st.title("Willkommen bei AI Primer Design Pro 🧬")
        st.markdown("""
        Willkommen bei **AI Primer Design Pro**,  
        deiner intelligenten Bioinformatik-Plattform für DNA-, RNA- und Protein-Analysen.  
        Hier kombinieren sich **KI**, **Laborautomatisierung** und **visuelle Werkzeuge**,  
        um Forschungsprozesse zu vereinfachen und zu beschleunigen.
        """)
        st.info("🌗 Tipp: Du kannst im Seitenmenü zwischen **Dark- und Light-Mode** wechseln.")
    else:
        st.title("Welcome to AI Primer Design Pro 🧬")
        st.markdown("""
        Welcome to **AI Primer Design Pro**,  
        your intelligent bioinformatics platform for DNA, RNA, and protein analysis.  
        Combining **AI**, **automation**, and **visual lab tools**  
        to simplify and accelerate research workflows.
        """)
        st.info("🌗 Tip: You can switch between **Dark and Light mode** in the sidebar.")

# --- MODULE: Sequence Management ---
elif menu == "🧬 Sequence Management":
    run_sequence_management()

# --- MODULE: Primer Design ---
elif menu == "🧫 Primer Design":
    run_primer_design()

# --- MODULE: Primer Design – Advanced ---
elif menu == "🧪 Primer Design – Advanced":
    run_primer_design_advanced()

# --- MODULE: Visual DNA Map ---
elif menu == "🧬 Visual DNA Map":
    import plotly.graph_objects as go
    from Bio import SeqIO
    import numpy as np
    import random

    st.title("🧬 Geneious-Style Visual DNA Map")

    # UI-Optionen
    colA, colB, colC = st.columns(3)
    with colA:
        show_features = st.checkbox("🧩 Features (Promoter/CDS) anzeigen", value=True)
    with colB:
        show_heatmap = st.checkbox("🌡️ Heatmap (GC% / AI-Score) anzeigen", value=True)
    with colC:
        detect_features = st.checkbox("🧠 Automatische Feature-Erkennung (GenBank)", value=True)

    uploaded_file = st.file_uploader("📂 Datei hochladen (FASTA oder GenBank)", type=["fasta", "fa", "gb", "gbk"])

    if not uploaded_file:
        st.info("⬆️ Bitte lade eine Sequenzdatei hoch (.fasta oder .gbk), um die DNA-Karte zu generieren.")
    else:
        # Datei speichern
        file_path = "uploaded_sequence.tmp"
        with open(file_path, "wb") as f:
            f.write(uploaded_file.getbuffer())

        # GenBank oder FASTA einlesen
        record = None
        if uploaded_file.name.endswith((".gb", ".gbk")):
            from Bio import SeqIO
            record = SeqIO.read(file_path, "genbank")
        else:
            record = SeqIO.read(file_path, "fasta")
        seq_length = len(record.seq)

        # --- Automatische Feature-Erkennung aus GenBank ---
        features = []
        if detect_features and hasattr(record, "features"):
            color_map = {
                "gene": "#1f77b4",
                "CDS": "#2ca02c",
                "promoter": "#ff7f0e",
                "misc_feature": "#9467bd",
                "tRNA": "#8c564b",
                "rRNA": "#e377c2",
                "exon": "#17becf"
            }
            for f in record.features:
                ftype = f.type.lower()
                if ftype in color_map and "location" in dir(f):
                    start = int(f.location.start)
                    end = int(f.location.end)
                    features.append({
                        "name": ftype.upper(),
                        "start": start,
                        "end": end,
                        "color": color_map[ftype]
                    })
        else:
            # Default Beispiel-Daten
            features = [
                {"name": "Promoter", "start": 20, "end": 60, "color": "#ffa600"},
                {"name": "CDS", "start": 80, "end": 480, "color": "#66b3ff"},
            ]

        # Beispiel-Primer (später KI-generiert)
        primers = [
            {"name": "Fwd1", "start": 120, "end": 140, "Tm": 59.2, "GC": 45, "score": 90, "strand": "+"},
            {"name": "Rev1", "start": 420, "end": 440, "Tm": 61.5, "GC": 52, "score": 72, "strand": "-"},
        ]

        # --- Plot ---
        fig = go.Figure()

        # 1️⃣ Heatmap (schmal unter der DNA)
        if show_heatmap:
            heat = np.array([random.randint(40, 70) for _ in range(seq_length)], dtype=float)
            z = np.expand_dims(heat, axis=0)
            fig.add_trace(go.Heatmap(
                z=z,
                x=list(range(seq_length)),
                y=[-0.25],
                colorscale="Viridis",
                opacity=0.35,
                showscale=True,
                name="GC% / AI Heatmap",
                hoverinfo="x+z",
                colorbar=dict(thickness=12, len=0.5, y=0.7)
            ))

        # 2️⃣ DNA-Basislinie
        fig.add_trace(go.Scatter(
            x=[0, seq_length],
            y=[0, 0],
            mode="lines",
            line=dict(color="#bfbfbf", width=12),
            name="DNA",
            hoverinfo="skip"
        ))

        # 3️⃣ Features
        if show_features and features:
            for f in features:
                fig.add_trace(go.Scatter(
                    x=[f["start"], f["end"]],
                    y=[0, 0],
                    mode="lines",
                    line=dict(color=f["color"], width=16),
                    name=f["name"],
                    hovertemplate=f"{f['name']}<br>{f['start']}–{f['end']} bp<extra></extra>"
                ))

        # 4️⃣ Primer
        for p in primers:
            color = "#00cc00" if p["strand"] == "+" else "#ff4d4d"
            arrow = "▶" if p["strand"] == "+" else "◀"
            fig.add_trace(go.Scatter(
                x=[p["start"], p["end"]],
                y=[0, 0],
                mode="lines+text",
                line=dict(color=color, width=20, shape="hv"),
                text=f"{arrow} {p['name']}",
                textposition="top center",
                name=p["name"],
                hovertemplate=(
                    f"<b>{p['name']}</b><br>"
                    f"Pos {p['start']}–{p['end']} bp<br>"
                    f"Tm {p['Tm']} °C · GC {p['GC']} %<br>"
                    f"AI-Score {p['score']}<extra></extra>"
                )
            ))

        # Layout
        fig.update_layout(
            title="🧬 Geneious-Style Visual DNA Map with Primer Heatmap & Auto-Features",
            xaxis_title="Nukleotidposition (bp)",
            yaxis=dict(visible=False, range=[-1, 1]),
            showlegend=True,
            height=420,
            plot_bgcolor="white",
            margin=dict(l=20, r=20, t=60, b=20),
            hovermode="x unified"
        )

        st.success("✅ DNA-Karte generiert!")
        st.plotly_chart(fig, use_container_width=True)

# --- MODULE: Cloning & Assembly Tools ---
elif menu == "🧫 Cloning & Assembly Tools":
    run_cloning_tools()

# --- MODULE: Protein Tools ---
elif menu == "🧬 Protein Tools":
    run_protein_tools()

# --- MODULE: Plasmid Designer ---
elif menu == "🧫 Plasmid Designer":
    run_plasmid_designer()

# --- MODULE: Plasmid Plus ---
elif menu == "🧬 Plasmid Plus":
    run_plasmid_plus()

# --- MODULE: Database & Reference Integration ---
elif menu == "🧫 Database & Reference Integration":
    run_database_integration()

# --- MODULE: Data Management ---
elif menu == "🧬 Data Management":
    run_data_management()

# --- MODULE: Alignment & Phylogeny ---
elif menu == "🌳 Alignment & Phylogeny":
    run_alignment_phylogeny()

# --- MODULE: AI Learning & Chatbot System ---
elif menu == "🤖 AI Learning & Chatbot System":
    run_ai_learning_chatbot()

# --- MODULE: Reports & Export Center ---
elif menu == "📊 Reports & Export Center":
    run_reports_export_center()

# --- Footer ---
st.markdown("---")
st.caption("🧠 Entwickelt mit ❤️ in Hamburg · Version 3.0 · Visual DNA Map + Auto Feature Detection · Zweisprachig DE/EN")
