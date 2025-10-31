# -*- coding: utf-8 -*-
"""
Geneious-Style Visual DNA Map (v3.2)
Erweiterte Darstellung mit:
- FASTA- oder GenBank-Upload ODER Direkteingabe
- GC%-Heatmap (Fensterweise berechnet)
- Feature-Anzeige (Promoter, CDS)
- Primer-Farben (Fwd/Rev)
"""

import streamlit as st
import plotly.graph_objects as go
from Bio import SeqIO
from io import StringIO

# =====================================================
# 🔬 Hilfsfunktionen
# =====================================================

def calc_gc_windows(sequence, window=10):
    """Berechne GC% pro Fenster für Heatmap."""
    seq = sequence.upper()
    gc_values = []
    for i in range(0, len(seq), window):
        win = seq[i:i + window]
        gc = (win.count("G") + win.count("C")) / len(win) * 100 if len(win) > 0 else 0
        gc_values.append(gc)
    return gc_values


def draw_visual_dna_map(seq_record, primers=None, show_heatmap=True, show_features=True):
    """
    Erstellt interaktive DNA-Karte mit Primern, GC%-Heatmap und Features.
    """
    seq = str(seq_record.seq)
    seq_len = len(seq)

    # ---------- GC%-Heatmap ----------
    x_heat = []
    y_heat = []
    color_heat = []

    if show_heatmap:
        gc_values = calc_gc_windows(seq)
        for i, gc in enumerate(gc_values):
            start = i * 10
            end = start + 10
            x_heat.extend([start, end, None])
            y_heat.extend([0, 0, None])
            # Farbverlauf (mehr GC = dunkler)
            color_heat.append(f"rgba({255 - int(gc*2)}, {255 - int(gc)}, {180}, 0.35)")

    # ---------- Figure Basis ----------
    fig = go.Figure()

    # DNA-Basislinie
    fig.add_trace(go.Scatter(
        x=[0, seq_len],
        y=[0, 0],
        mode="lines",
        line=dict(color="lightgray", width=8),
        name="DNA"
    ))

    # Heatmap Layer (schmale Bänder)
    if show_heatmap:
        fig.add_trace(go.Scatter(
            x=[i * 10 for i in range(len(calc_gc_windows(seq)))],
            y=[-0.05] * len(calc_gc_windows(seq)),
            mode="markers",
            marker=dict(
                size=12,
                color=calc_gc_windows(seq),
                colorscale="YlOrBr",
                showscale=True,
                colorbar=dict(title="GC%", thickness=10)
            ),
            name="GC%-Heatmap"
        ))

    # ---------- Primer Layer ----------
    if primers:
        for p in primers:
            color = "green" if "Fwd" in p["name"] else "red"
            fig.add_trace(go.Scatter(
                x=[p["start"], p["end"]],
                y=[0, 0],
                mode="lines",
                line=dict(width=12, color=color),
                name=p["name"],
                hovertext=(
                    f"<b>{p['name']}</b><br>"
                    f"Position: {p['start']}–{p['end']}<br>"
                    f"Tm: {p.get('Tm', '–')}°C<br>"
                    f"GC: {p.get('GC', '–')}%<br>"
                    f"Score: {p.get('score', '–')}"
                ),
                hoverinfo="text"
            ))

    # ---------- Features aus GenBank ----------
    if show_features and hasattr(seq_record, "features"):
        for f in seq_record.features:
            try:
                start = int(f.location.start)
                end = int(f.location.end)
                feature_name = f.qualifiers.get("gene", ["Feature"])[0]
                ftype = f.type.upper()
                color = "royalblue" if "PROMOTER" in ftype else "orange"

                fig.add_trace(go.Scatter(
                    x=[start, end],
                    y=[0.15, 0.15],
                    mode="lines",
                    line=dict(width=8, color=color),
                    name=f"{ftype}: {feature_name}",
                    hovertext=f"{ftype} ({feature_name})<br>{start}-{end} bp",
                    hoverinfo="text"
                ))
            except Exception:
                pass

    # ---------- Layout ----------
    fig.update_layout(
        title="🧬 Geneious-Style Visual DNA Map (thin-band Heatmap + Auto-Features)",
        showlegend=True,
        xaxis_title="Position (bp)",
        yaxis_visible=False,
        plot_bgcolor="white",
        height=400,
        margin=dict(l=40, r=40, t=60, b=40)
    )

    st.success("✅ DNA-Karte generiert!")
    st.plotly_chart(fig, use_container_width=True)


# =====================================================
# 🚀 Streamlit App
# =====================================================

def visualize_dna_map():
    st.title("🧬 Geneious-Style Visual DNA Map")
    st.caption("Visualisiert DNA-Sequenzen mit GC%-Heatmap, Features und Primern")

    # Auswahl: Upload oder Eingabe
    mode = st.radio("🔹 Eingabeart wählen:", ["Datei hochladen", "Sequenz manuell eingeben"])

    seq_record = None
    if mode == "Datei hochladen":
        upl = st.file_uploader("FASTA oder GenBank-Datei hochladen", type=["fasta", "fa", "gb", "gbk"])
        if upl:
            try:
                if upl.name.lower().endswith((".gb", ".gbk")):
                    seq_record = SeqIO.read(upl, "genbank")
                else:
                    seq_record = SeqIO.read(upl, "fasta")
                st.success(f"✅ Datei erfolgreich geladen: {upl.name}")
            except Exception as e:
                st.error(f"Fehler beim Einlesen: {e}")

    else:
        seq_text = st.text_area("DNA-Sequenz (5′→3′, ACGT)", height=150)
        if st.button("✅ Sequenz übernehmen"):
            from Bio.Seq import Seq
            from Bio.SeqRecord import SeqRecord
            seq_record = SeqRecord(Seq(seq_text.strip()), id="ManualInput")

    if seq_record:
        primers = [
            {"name": "Fwd1", "start": 10, "end": 30, "Tm": 59.2, "GC": 45, "score": 88},
            {"name": "Rev1", "start": len(seq_record.seq) - 40, "end": len(seq_record.seq) - 20, "Tm": 60.5, "GC": 52, "score": 76},
        ]
        draw_visual_dna_map(seq_record, primers, show_heatmap=True, show_features=True)
