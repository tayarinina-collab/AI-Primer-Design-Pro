import streamlit as st
import plotly.graph_objects as go

def render():
    st.header("Plasmid-Karte (Basis)")
    st.caption("Basic circular map preview / Einfache zirkulaere Karten-Vorschau")

    plasmid_len = st.number_input("Plasmid length (bp)", 500, 20000, 3000, step=100)
    features = st.text_area("Features (CSV: name,start,end)", height=120,
                            value="ori,100,300\nAmpR,800,1600\nGFP,1800,2600")

    theta, r, text = [], [], []
    for line in features.splitlines():
        try:
            name, s, e = line.split(",")
            s = int(s); e = int(e)
            mid = (s+e)/2 / plasmid_len * 360.0
            theta.append(mid); r.append(1); text.append(f"{name} ({s}-{e})")
        except Exception:
            continue

    fig = go.Figure(go.Scatterpolar(theta=theta, r=r, mode="markers+text", text=text, textposition="top center"))
    fig.update_layout(polar=dict(radialaxis=dict(visible=False), angularaxis=dict(dtick=45)),
                      showlegend=False, margin=dict(l=10,r=10,t=10,b=10))

    st.plotly_chart(fig, use_container_width=True)
    st.info("This is a basic visualization. Advanced maps (enzymes, annotations) planned via BioPython Graphics.")