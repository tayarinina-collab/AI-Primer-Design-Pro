# Title: 🧬 Plasmid Plus
# -*- coding: utf-8 -*-
"""
Plasmid Plus – Erweiterte Ansicht
- Analyse von Feature-Längen
- Prozentuale Feature-Berechnung
- GenBank + CSV Export
- Vorschau mit Farbe und Statistik
"""

import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
from io import StringIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
import math


def _parse_features(text):
    rows = []
    for ln in text.splitlines():
        ln = ln.strip()
        if not ln or ln.startswith("#"):
            continue
        parts = [p.strip() for p in ln.split(",")]
        if len(parts) < 3:
            continue
        name, s, e = parts[:3]
        try:
            s, e = int(s), int(e)
            rows.append({"name": name, "start": s, "end": e, "length": e - s})
        except:
            pass
    return pd.DataFrame(rows)


def _plot_features(df, plasmid_len):
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(6, 6))
    ax.set_theta_direction(-1)
    ax.set_theta_offset(math.pi / 2)
    ax.set_xticks([])
    ax.set_yticks([])

    # Circle
    theta = [math.radians(t) for t in range(360)]
    ax.plot(theta, [1]*360, color="gray", alpha=0.4)

    # Features
    colors = ["tab:blue", "tab:green", "tab:red", "tab:orange", "tab:purple"]
    for i, r in df.iterrows():
        start = 2 * math.pi * r["start"] / plasmid_len
        end = 2 * math.pi * r["end"] / plasmid_len
        if end < start:
            end += 2 * math.pi
        color = colors[i % len(colors)]
        ts = [start + (end - start) * k / 100 for k in range(101)]
        ax.plot(ts, [0.9]*len(ts), linewidth=8, color=color, alpha=0.4)
        ax.text((start+end)/2, 1.05, r["name"], ha="center", va="center", color=color, fontsize=9)
    st.pyplot(fig, use_container_width=True)


def run_plasmid_plus():
    st.title("🧬 Plasmid Plus – Erweiterte Analyse")

    c1, c2 = st.columns(2)
    with c1:
        plasmid_len = st.number_input("🔢 Plasmidlänge (bp)", 500, 50000, 5000, step=100)
    with c2:
        st.markdown(" ")

    st.markdown("**Features (CSV: name,start,end)**")
    example = "ori,100,300\nAmpR,700,1600\nGFP,1800,2600\nMCS,3000,3300"
    features = st.text_area("Feature-Liste", value=example, height=120)
    df = _parse_features(features)

    if not df.empty:
        st.dataframe(df, use_container_width=True, hide_index=True)

        # Statistik
        total_bp = df["length"].sum()
        percent = round(100 * total_bp / plasmid_len, 2)
        st.markdown(f"**Gesamt Feature-Länge:** {total_bp} bp ({percent}% des Plasmids)")

        # Plot
        _plot_features(df, plasmid_len)

        # Export
        seq = Seq("ATGC" * (plasmid_len // 4))
        rec = SeqRecord(seq, id="plasmid_plus", description="AI Primer Design Pro – Plasmid Plus")
        for _, r in df.iterrows():
            rec.features.append(
                SeqFeature(
                    FeatureLocation(r["start"], r["end"]),
                    type="misc_feature",
                    qualifiers={"label": r["name"]}
                )
            )

        buf = StringIO()
        SeqIO.write(rec, buf, "genbank")
        gb_bytes = buf.getvalue().encode("utf-8")

        st.download_button("💾 GenBank (.gb) exportieren", data=gb_bytes,
                           file_name="plasmid_plus.gb", mime="application/octet-stream")
        st.download_button("📄 CSV exportieren", data=df.to_csv(index=False).encode("utf-8"),
                           file_name="plasmid_plus_features.csv", mime="text/csv")
    else:
        st.info("Bitte Features eingeben, um Karte und Analyse zu sehen.")
