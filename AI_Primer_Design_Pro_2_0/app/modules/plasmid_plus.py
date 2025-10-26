# -*- coding: utf-8 -*-
"""
Plasmid Designer – Basis + GenBank-Export
- Kachel: einfache kreisförmige Map
- Feature-Editor (CSV)
- SeqRecord + GenBank-Export (.gb)
- CSV-Export der Features
"""

import io
import math
import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation


def _parse_features_csv(text: str) -> pd.DataFrame:
    """
    Erwartet CSV-Zeilen: name,start,end[,type[,strand]]
    - start/end: 1-basiert (wie üblich in GenBank). Intern wird 0-basiert verwendet.
    - type (optional): default 'misc_feature'
    - strand (optional): +1, -1 oder leer => None
    """
    rows = []
    for ln in text.splitlines():
        ln = ln.strip()
        if not ln or ln.startswith("#"):
            continue
        parts = [p.strip() for p in ln.split(",")]
        if len(parts) < 3:
            continue
        name, s, e = parts[:3]
        ftype = parts[3] if len(parts) >= 4 and parts[3] else "misc_feature"
        strand = parts[4] if len(parts) >= 5 and parts[4] else ""
        try:
            s = int(s)  # 1-based
            e = int(e)
            strand = int(strand) if strand in ("1", "-1") else 0
            rows.append(
                {"name": name, "start_1b": s, "end_1b": e, "type": ftype, "strand": strand}
            )
        except Exception:
            # ignorieren fehlerhafte Zeilen
            pass
    if not rows:
        return pd.DataFrame(columns=["name", "start_1b", "end_1b", "type", "strand"])
    df = pd.DataFrame(rows)
    return df


def _make_seqrecord(plasmid_len: int, features_df: pd.DataFrame, circular=True) -> SeqRecord:
    """
    Erzeugt einen SeqRecord mit Dummy-DNA-Sequenz (A/C/G/T), Features & Qualifiers.
    Reihenfolge/Anzahl der Basen ist hier egal – die Features/Annotationen sind entscheidend
    für Tools wie Geneious/SnapGene.
    """
    # Dummy-DNA (reproduzierbar): wiederholendes Muster
    bases = ("ATGCGT" * ((plasmid_len // 6) + 1))[:plasmid_len]
    seq = Seq(bases)

    rec = SeqRecord(
        seq=seq,
        id="plasmid_1",
        name="Plasmid_1",
        description="AI Primer Design Pro – Synthetic plasmid",
    )
    # Metadata für Kreisförmigkeit
    rec.annotations["molecule_type"] = "DNA"
    rec.annotations["topology"] = "circular" if circular else "linear"
    rec.annotations["data_file_division"] = "SYN"
    rec.annotations["date"] = ""

    # Features (0-basiert intern!)
    for _, r in features_df.iterrows():
        start0 = max(0, int(r["start_1b"]) - 1)
        end0 = min(plasmid_len, int(r["end_1b"]))  # GenBank-Endposition ist exklusiv in FeatureLocation
        ftype = str(r.get("type", "misc_feature"))
        strand = int(r.get("strand", 0))
        qual = {"label": str(r["name"]), "note": [str(r["name"])]}

        feat = SeqFeature(
            FeatureLocation(start0, end0, strand=strand if strand in (-1, 1) else 0),
            type=ftype,
            qualifiers=qual,
        )
        rec.features.append(feat)

    return rec


def _plot_circular_map(plasmid_len: int, features_df: pd.DataFrame):
    """Sehr einfache Polarplot-Vorschau (kein SnapGene-Ersatz, aber hilfreich)."""
    if plasmid_len <= 0:
        st.warning("Plasmidlänge muss > 0 sein.")
        return

    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(6, 6))
    ax.set_theta_direction(-1)
    ax.set_theta_offset(math.pi / 2.0)
    ax.set_xticks([])
    ax.set_yticks([])

    # Circle
    thetas = [t * math.pi / 180 for t in range(0, 360)]
    rs = [1] * len(thetas)
    ax.plot(thetas, rs, linewidth=2, alpha=0.4)

    # Features als Bögen
    colors = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple", "tab:brown"]
    for i, r in features_df.iterrows():
        s = int(r["start_1b"])  # 1-based
        e = int(r["end_1b"])
        if s < 1 or e < 1 or s > plasmid_len or e > plasmid_len:
            continue
        # Start/Ende -> Winkel
        theta_s = 2 * math.pi * (s - 1) / plasmid_len
        theta_e = 2 * math.pi * (e - 1) / plasmid_len
        col = colors[i % len(colors)]
        # Arc zeichnen (r leicht <1, damit es sichtbar ist)
        steps = 120
        if theta_e < theta_s:
            theta_e += 2 * math.pi
        ts = [theta_s + (theta_e - theta_s) * k / steps for k in range(steps + 1)]
        ax.plot(ts, [0.92] * (steps + 1), color=col, linewidth=8, alpha=0.35)
        # Label
        mid = (theta_s + theta_e) / 2
        ax.text(mid, 1.06, str(r["name"]), ha="center", va="center", fontsize=9, color=col)

    st.pyplot(fig, use_container_width=True)


def run_plasmid_designer():
    st.title("🧫 Plasmid-Karte (Basis) + GenBank-Export")

    c1, c2 = st.columns([1, 1])
    with c1:
        plasmid_len = st.number_input("Plasmidlänge (bp)", 100, 200000, 5000, step=100)
    with c2:
        circular = st.toggle("Topologie: Kreisförmig", value=True)

    st.markdown("**Features (CSV: name,start,end[,type[,strand]])**  – Beispiele:")
    with st.expander("📋 Beispiel einblenden", expanded=False):
        st.code("ori,100,300\nAmpR,800,1600,gene,1\nGFP,1800,2600,CDS,1\nMCS,2650,2750,misc_feature,0", language="text")

    features_csv = st.text_area(
        "Feature-Liste",
        value="ori,100,300\nAmpR,800,1600\ngfp,1800,2600,CDS,1\nMCS,2650,2750",
        height=140,
    )

    # Parsen + Anzeigen
    df = _parse_features_csv(features_csv)
    st.dataframe(df, use_container_width=True, hide_index=True)

    st.markdown("---")
    st.subheader("🗺️ Vorschau (schematisch)")
    _plot_circular_map(plasmid_len, df)

    st.markdown("---")
    st.subheader("⬇️ Export")

    # GenBank-Export
    rec = _make_seqrecord(plasmid_len, df, circular=circular)
    gb_buf = io.StringIO()
    from Bio import SeqIO
    SeqIO.write(rec, gb_buf, "genbank")
    gb_bytes = gb_buf.getvalue().encode("utf-8")

    st.download_button(
        "💾 GenBank (.gb) exportieren",
        data=gb_bytes,
        file_name="plasmid.gb",
        mime="application/octet-stream",
    )

    # CSV-Export der Features
    st.download_button(
        "📄 Features als CSV exportieren",
        data=df.to_csv(index=False).encode("utf-8"),
        file_name="plasmid_features.csv",
        mime="text/csv",
    )

    st.caption("Hinweis: Die GenBank-Datei enthält die Features/Annotationen und eine synthetische DNA-Sequenz.")
