import streamlit as st
import plotly.graph_objects as go
from Bio import SeqIO
from io import StringIO

def visualize_dna_map():
    """
    Interaktive DNA-Map (Geneious-Stil)
    - Eingabe: FASTA-Datei oder direkte Sequenzeingabe
    - Darstellung: DNA-Linie + Primerpositionen
    """

    st.title("🧬 Visual DNA Map – Geneious Style")
    st.caption("Interaktive Darstellung der DNA-Sequenz mit markierten Primern")

    # ================= Eingabeoptionen =================
    st.subheader("📥 Eingabeoptionen")
    col1, col2 = st.columns(2)
    with col1:
        uploaded_file = st.file_uploader("FASTA-Datei hochladen", type=["fasta", "fa", "txt"])
    with col2:
        seq_input = st.text_area("Oder Sequenz direkt eingeben (5'→3')", height=150,
                                 placeholder=">optional_header\nATGAAAGCTGTTGCTGCTG...")

    # ================= Sequenz einlesen =================
    seq = ""
    if uploaded_file is not None:
        try:
            content = uploaded_file.getvalue().decode("utf-8")
            record = SeqIO.read(StringIO(content), "fasta")
            seq = str(record.seq).upper()
        except Exception as e:
            st.error(f"❌ Fehler beim Lesen der FASTA-Datei: {e}")
            return
    elif seq_input.strip():
        txt = seq_input.strip()
        if txt.startswith(">"):
            # FASTA-Format in Textbox erkannt
            try:
                record = SeqIO.read(StringIO(txt), "fasta")
                seq = str(record.seq).upper()
            except Exception:
                seq = "".join([c for c in txt if c in "ACGT"])
        else:
            seq = "".join([c for c in txt.upper() if c in "ACGT"])
    else:
        st.info("Bitte lade eine FASTA-Datei hoch oder gib eine DNA-Sequenz ein.")
        return

    seq_length = len(seq)
    st.success(f"✅ Sequenz geladen ({seq_length} bp)")

    # ================= Dummy-Primer (Beispiel) =================
    primers = [
        {"name": "Fwd1", "start": 30, "end": 50, "Tm": 59.2, "GC": 45, "score": 88},
        {"name": "Rev1", "start": 120, "end": 140, "Tm": 60.5, "GC": 52, "score": 72},
        {"name": "Fwd2", "start": 200, "end": 220, "Tm": 61.3, "GC": 50, "score": 90},
        {"name": "Rev2", "start": 300, "end": 320, "Tm": 58.9, "GC": 47, "score": 80},
    ]

    # ================= Plotly Visualisierung =================
    fig = go.Figure()

    # DNA-Strang
    fig.add_trace(go.Scatter(
        x=[0, seq_length],
        y=[0, 0],
        mode="lines",
        line=dict(color="lightgray", width=6),
        name="DNA"
    ))

    # Primer hinzufügen
    for p in primers:
        color = f"rgba({255 - int(p['score']*2.55)}, {int(p['score']*2.55)}, 0, 0.8)"
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
                f"Score: {p.get('score', '–')}"
            ),
            hoverinfo="text"
        ))

    # Layout-Einstellungen
    fig.update_layout(
        title="🧬 Visual DNA Map mit Primer-Markierung",
        showlegend=True,
        xaxis_title="Position (bp)",
        yaxis_visible=False,
        plot_bgcolor="white",
        height=400,
        margin=dict(l=20, r=20, t=50, b=20)
    )

    st.plotly_chart(fig, use_container_width=True)

    # ================= Download =================
    st.download_button(
        "⬇️ Sequenz als FASTA exportieren",
        data=f">visual_dna_map_sequence\n{seq}\n",
        file_name="visual_dna_map_sequence.fasta"
    )
