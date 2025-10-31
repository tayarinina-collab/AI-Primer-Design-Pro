import streamlit as st
import plotly.graph_objects as go
from Bio import SeqIO

def visualize_dna_map(fasta_path, primers, color_by="type"):
    """
    Zeigt eine interaktive DNA-Karte mit markierten Primern.
    primers: Liste mit Dictionaries, z. B.:
    [
        {"name": "Fwd1", "start": 120, "end": 140, "Tm": 59.2, "GC": 45, "score": 88},
        {"name": "Rev1", "start": 420, "end": 440, "Tm": 60.5, "GC": 52, "score": 72}
    ]
    """

    # FASTA-Sequenz lesen
    record = SeqIO.read(fasta_path, "fasta")
    seq_length = len(record.seq)

    # Baseline für DNA
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=[0, seq_length],
        y=[0, 0],
        mode="lines",
        line=dict(color="lightgray", width=6),
        name="DNA"
    ))

    # Primer als Balken mit Farbe/Heatmap
    for p in primers:
        if color_by == "score":
            # AI-Score in Farbe umwandeln (0–100)
            color = f"rgba({255 - int(p['score']*2.55)}, {int(p['score']*2.55)}, 0, 0.8)"
        else:
            color = "red" if "Fwd" in p["name"] else "blue"

        fig.add_trace(go.Scatter(
            x=[p["start"], p["end"]],
            y=[0, 0],
            mode="lines",
            line=dict(width=10, color=color),
            name=p["name"],
            hovertext=(
                f"<b>{p['name']}</b><br>"
                f"Position: {p['start']}–{p['end']}<br>"
                f"Tm: {p['Tm']}°C<br>"
                f"GC: {p['GC']}%<br>"
                f"AI-Score: {p.get('score','–')}"
            )
        ))

    fig.update_layout(
        title="🧬 Visual DNA Map with Primer Heatmap",
        showlegend=True,
        xaxis_title="Nukleotidposition",
        yaxis_visible=False,
        plot_bgcolor="white",
        height=350
    )

    st.plotly_chart(fig, use_container_width=True)
