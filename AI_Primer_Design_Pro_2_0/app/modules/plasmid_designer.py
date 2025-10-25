import streamlit as st
import plotly.graph_objects as go

def run_plasmid_designer():
    st.header("🧬 Plasmid-Karte (Basis)")
    st.caption("Basic circular map preview / Einfache zirkuläre Karten-Vorschau")

    # Eingabe: Plasmidgröße
    plasmid_len = st.number_input("Plasmidlänge (bp)", 500, 20000, 3000, step=100)

    # Eingabe: Gen-/Feature-Positionen
    features = st.text_area(
        "Features (CSV: Name, Start, Ende)",
        height=120,
        value="ori,100,300\nAmpR,800,1600\nGFP,1800,2600"
    )

    # Daten vorbereiten
    theta, r, text = [], [], []
    for line in features.splitlines():
        try:
            name, s, e = line.split(",")
            s = int(s)
            e = int(e)
            mid = (s + e) / 2 / plasmid_len * 360.0
            theta.append(mid)
            r.append(1)
            text.append(f"{name} ({s}-{e})")
        except Exception:
            continue

    # Plotly-Kreisdiagramm
    fig = go.Figure(
        go.Scatterpolar(
            theta=theta,
            r=r,
            mode="markers+text",
            text=text,
            textposition="top center",
            marker=dict(size=12, color="lightblue", line=dict(color="darkblue", width=2))
        )
    )

    # Layout anpassen
    fig.update_layout(
        polar=dict(
            radialaxis=dict(visible=False),
            angularaxis=dict(dtick=45)
        ),
        showlegend=False,
        margin=dict(l=10, r=10, t=10, b=10),
        template="plotly_dark" if st.session_state.get("theme", "light") == "dark" else "plotly_white"
    )

    # Anzeige
    st.plotly_chart(fig, use_container_width=True)
    st.info("Dies ist eine Basis-Visualisierung. Erweiterte Karten (Enzyme, Annotationen) werden über BioPython-Integration ergänzt.")
