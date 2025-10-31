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
        Willkommen bei **AI Primer Design Pro** – deiner intelligenten Bioinformatik-Plattform
        für DNA-, RNA- und Protein-Analysen. Kombiniert **KI**, **Automatisierung** und **visuelle Tools**,
        um Forschungsprozesse zu beschleunigen.
        """)
        st.info("🌗 Tipp: Im Seitenmenü zwischen **Dark- und Light-Mode** umschalten.")
    else:
        st.title("Welcome to AI Primer Design Pro 🧬")
        st.markdown("""
        Welcome to **AI Primer Design Pro** – your intelligent bioinformatics platform for DNA, RNA,
        and protein analysis. It combines **AI**, **automation**, and **visual tools** to accelerate research.
        """)
        st.info("🌗 Tip: Switch **Dark/Light mode** in the sidebar.")

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
        show_heatmap = st.checkbox("🌡️ Heatmap (GC%/AI) anzeigen", value=True)
    with colC:
        detect_features = st.checkbox("🧠 Auto-Features aus GenBank", value=True)

    uploaded_file = st.file_uploader("📂 Datei hochladen (FASTA/GenBank)", type=["fasta", "fa", "gb", "gbk"])

    if not uploaded_file:
        st.info("⬆️ Bitte lade eine Sequenzdatei (.fasta / .gbk) hoch.")
    else:
        # Datei speichern & laden
        file_path = "uploaded_sequence.tmp"
        with open(file_path, "wb") as f:
            f.write(uploaded_file.getbuffer())

        # Record parsen
        if uploaded_file.name.endswith((".gb", ".gbk")):
            record = SeqIO.read(file_path, "genbank")
        else:
            record = SeqIO.read(file_path, "fasta")

        seq_length = len(record.seq)

        # ---------- Helper: kompakter GC%-Score (Sliding Window) ----------
        def sliding_gc_scores(seq, window=15):
            seq = str(seq).upper()
            if window < 2:
                window = 2
            scores = []
            half = window // 2
            for i in range(seq_length):
                a = max(0, i - half)
                b = min(seq_length, i + half)
                sub = seq[a:b]
                gc = (sub.count("G") + sub.count("C")) / max(1, len(sub))
                scores.append(gc)  # 0..1
            return np.array(scores, dtype=float)

        # ---------- Heatmap-Farben (schmale Shapes unter DNA) ----------
        def color_from_score(s):
            # einfache Viridis-nahe Mischung: 0→lila, 1→gelbgrün
            # Alpha 0.35 für Hintergrund
            r = int(68 + s * (253 - 68))   # grobe Approx an Viridis
            g = int(1 + s * (231 - 1))
            b = int(84 + s * (37 - 84))
            return f"rgba({r},{g},{b},0.35)"

        # ---------- Features ----------
        features = []
        if detect_features and hasattr(record, "features"):
            color_map = {
                "gene": "#1f77b4",
                "cds": "#2ca02c",
                "promoter": "#ff7f0e",
                "misc_feature": "#9467bd",
                "trna": "#8c564b",
                "rrna": "#e377c2",
                "exon": "#17becf",
            }
            for feat in getattr(record, "features", []):
                ftype = feat.type.lower()
                if ftype in color_map:
                    start = int(feat.location.start)
                    end = int(feat.location.end)
                    if start < end:
                        features.append({"name": ftype.upper(), "start": start, "end": end, "color": color_map[ftype]})
        else:
            # Fallback-Beispiele
            features = [
                {"name": "Promoter", "start": 20, "end": 60, "color": "#ffb000"},
                {"name": "CDS", "start": 80, "end": min(480, seq_length-1), "color": "#66b3ff"},
            ]

        # ---------- Demo-Primer (später aus Design übernehmen) ----------
        primers = [
            {"name": "Fwd1", "start": max(0, int(seq_length*0.25)-10), "end": max(1, int(seq_length*0.25)+10),
             "Tm": 59.2, "GC": 45, "score": 90, "strand": "+"},
            {"name": "Rev1", "start": max(0, int(seq_length*0.9)-10), "end": max(1, int(seq_length*0.9)+10),
             "Tm": 61.5, "GC": 52, "score": 72, "strand": "-"},
        ]

        # ---------- Plot ----------
        fig = go.Figure()

        # 1) Heatmap als SCHMALES Band via Shapes (kein go.Heatmap mehr)
        if show_heatmap:
            scores = sliding_gc_scores(record.seq, window=15)  # 0..1
            shapes = []
            y0, y1 = -0.25, 0.25  # dünnes Band
            for i, s in enumerate(scores):
                shapes.append(dict(
                    type="rect",
                    xref="x", yref="y",
                    x0=i, x1=i+1, y0=y0, y1=y1,
                    line=dict(width=0),
                    fillcolor=color_from_score(s),
                    layer="below",
                ))
            fig.update_layout(shapes=shapes)

        # 2) DNA-Linie
        fig.add_trace(go.Scatter(
            x=[0, seq_length], y=[0, 0],
            mode="lines", line=dict(color="#bfbfbf", width=12),
            name="DNA", hoverinfo="skip"
        ))

        # 3) Features (Promoter/CDS/…)
        if show_features and features:
            for f in features:
                fig.add_trace(go.Scatter(
                    x=[f["start"], f["end"]], y=[0, 0],
                    mode="lines",
                    line=dict(color=f["color"], width=16),
                    name=f["name"],
                    hovertemplate=f"{f['name']}<br>{f['start']}–{f['end']} bp<extra></extra>"
                ))

        # 4) Primer mit Richtung & Label
        for p in primers:
            color = "#00cc00" if p["strand"] == "+" else "#ff4d4d"
            arrow = "▶" if p["strand"] == "+" else "◀"
            fig.add_trace(go.Scatter(
                x=[p["start"], p["end"]], y=[0, 0],
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
            title="🧬 Geneious-Style Visual DNA Map (thin-band Heatmap + Auto-Features)",
            xaxis_title="Nukleotidposition (bp)",
            yaxis=dict(visible=False, range=[-1, 1]),
            showlegend=True,
            height=480,
            plot_bgcolor="white",
            margin=dict(l=20, r=20, t=60, b=20),
            hovermode="x unified",
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
st.caption("🧠 Entwickelt mit ❤️ in Hamburg · Version 3.1 · Visual DNA Map (thin-band) + Auto-Features")
