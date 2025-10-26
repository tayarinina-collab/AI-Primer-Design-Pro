# -*- coding: utf-8 -*-
"""
Cloning & Assembly Tools – extended
Gibson, Golden Gate, Restriktionsanalyse, Verdau, Plasmid-Map
+ Plasmid Annotator, Cloning-Primer Generator, Gateway/TA/TOPO, GenBank-Viewer
"""
from __future__ import annotations
import io
import re
from typing import List, Tuple, Dict

import streamlit as st
import pandas as pd
import plotly.graph_objects as go

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO, Restriction
from Bio.SeqFeature import SeqFeature, FeatureLocation

# optional primer3
try:
    import primer3
    P3_OK = True
except Exception:
    P3_OK = False

# ----------------------------- helpers ---------------------------------- #

def revcomp(s: str) -> str:
    return str(Seq(s).reverse_complement())

def gc(seq: str) -> float:
    s = seq.upper().replace("U", "T")
    return 0 if not s else 100.0 * (s.count("G") + s.count("C")) / len(s)

def linear_feature_plot(length: int, features: List[Tuple[int, int, str]]) -> go.Figure:
    """
    Simple linear map: features = [(start, end, label), ...]
    """
    fig = go.Figure()
    fig.add_shape(type="rect", x0=0, x1=length, y0=0.35, y1=0.65,
                  line=dict(color="#999"), fillcolor="#e8f5e9")
    for s, e, name in features:
        fig.add_shape(type="rect", x0=s, x1=e, y0=0.15, y1=0.85,
                      line=dict(color="#2e7d32"), fillcolor="#c8e6c9")
        fig.add_annotation(x=(s+e)/2, y=0.9, text=name, showarrow=False, font=dict(size=10))
    fig.update_xaxes(range=[0, length], title="Position (bp)")
    fig.update_yaxes(visible=False)
    fig.update_layout(height=200, margin=dict(l=40, r=20, t=10, b=30))
    return fig

def circular_map(length: int, features: List[Tuple[int, int, str]]) -> go.Figure:
    """
    Simple circular plasmid map using polar coordinates.
    """
    fig = go.Figure()
    # backbone
    fig.add_trace(go.Scatterpolar(r=[1, 1], theta=[0, 360], mode="lines",
                                  line=dict(width=8, color="#81c784"),
                                  hoverinfo="skip", showlegend=False))
    # features
    for s, e, name in features:
        theta0 = 360.0 * (s/length)
        theta1 = 360.0 * (e/length)
        fig.add_trace(go.Barpolar(r=[0.12], theta=[(theta0+theta1)/2],
                                  width=[max(2, theta1-theta0)],
                                  marker=dict(color="#558b2f"),
                                  name=name, hovertext=f"{name}: {s}..{e}", hoverinfo="text"))
    fig.update_layout(polar=dict(radialaxis=dict(visible=False),
                                 angularaxis=dict(direction='clockwise')),
                      showlegend=False, margin=dict(l=20, r=20, t=10, b=10),
                      height=350)
    return fig

# simple motif library for annotator
MOTIFS: Dict[str, List[str]] = {
    "ori": ["ori", "origin of replication", "colE1", "pMB1"],
    "AmpR": ["bla", "amp", "ampicillin", "beta-lactamase"],
    "KanR": ["kan", "kanamycin", "nptII", "neo"],
    "CatR": ["cat", "chloramphenicol"],
    "HisTag": ["histag", "his-tag"],
    "MCS": ["multiple cloning site", "mcs", "polylinker"],
    "LacZ": ["lacz", "alpha", "lacZα", "lacZa"],
    "T7 promoter": ["t7 promoter", "phi10"],
}

TYPE_IIS = {"BsaI": "GGTCTC", "BsmBI": "CGTCTC", "BbsI": "GAAGAC", "Esp3I": "CGTCTC"}

def search_motif(seq: str, pattern: str) -> List[Tuple[int, int]]:
    """naiver substring/regex matcher, gibt (start, end) zurück"""
    hits = []
    for m in re.finditer(pattern, seq, flags=re.IGNORECASE):
        hits.append((m.start(), m.end()))
    return hits

# --------------------------- main UI ------------------------------------ #

def run_cloning_tools():
    st.title("Cloning & Assembly Tools")
    st.caption("Gibson · Golden Gate · Restriktionen · Verdau/Ligation · Plasmid-Map · Annotator · Cloning-Primer · Gateway/TA/TOPO · GenBank Viewer")

    tabs = st.tabs([
        "Gibson",
        "Golden Gate",
        "Restriktionsanalyse",
        "Verdau",
        "Plasmid-Map",
        "Plasmid Annotator",
        "Cloning-Primer Generator",
        "Gateway / TA / TOPO",
        "GenBank-Viewer"
    ])

    # ---------- Gibson ----------
    with tabs[0]:
        st.subheader("Gibson Assembly Assistant")
        seq1 = st.text_area("Fragment 1 (5'→3')", height=120, key="gib1")
        seq2 = st.text_area("Fragment 2 (5'→3')", height=120, key="gib2")
        overlap = st.slider("Overlap (bp)", 10, 60, 20)
        if st.button("Gibson simulieren"):
            if seq1 and seq2:
                ov = seq1[-overlap:]
                if ov.upper() == seq2[:overlap].upper():
                    assembled = seq1 + seq2[overlap:]
                    st.success("Overlap erkannt – Contig erzeugt")
                    st.code(assembled)
                else:
                    st.warning("Kein passender Overlap gefunden.")
            else:
                st.info("Bitte beide Fragmente eingeben.")

    # ---------- Golden Gate ----------
    with tabs[1]:
        st.subheader("Golden Gate Designer")
        seq = st.text_area("DNA-Sequenz", height=120, key="gg_seq")
        enzyme = st.selectbox("Typ IIS Enzym", list(TYPE_IIS.keys()))
        if st.button("Schnittstellen suchen"):
            motif = TYPE_IIS[enzyme]
            pos = [i for i in range(len(seq) - len(motif) + 1) if seq[i:i+len(motif)].upper() == motif]
            if pos:
                st.success(f"{enzyme}: {len(pos)} Schnittstelle(n) bei {pos}")
            else:
                st.info("Keine Motive gefunden.")

    # ---------- Restriktionsanalyse ----------
    with tabs[2]:
        st.subheader("Restriktionsanalyse")
        up = st.file_uploader("FASTA laden", type=["fasta","fa","txt"], key="ra_up")
        if up:
            recs = list(SeqIO.parse(io.StringIO(up.getvalue().decode()), "fasta"))
            if not recs:
                st.warning("Keine Sequenz erkannt.")
            else:
                seq = recs[0].seq
                analysis = Restriction.Analysis(Restriction.AllEnzymes, seq)
                res = analysis.full()  # dict: {Enzym: [sites]}
                st.write(f"{len(res)} Enzyme mit Schnittstellen gefunden.")
                # kleine Tabelle
                rows = [{"Enzym": str(k), "Sites": ",".join(map(str, v))} for k, v in res.items() if v]
                st.dataframe(pd.DataFrame(rows), use_container_width=True)

    # ---------- Verdau ----------
    with tabs[3]:
        st.subheader("In-Silico Verdau")
        seq = st.text_area("DNA-Sequenz", height=120, key="vd_seq")
        enzyme = st.selectbox("Enzym", ["EcoRI","BamHI","HindIII","NotI"], key="vd_enzyme")
        if st.button("Verdau simulieren"):
            try:
                cut = getattr(Restriction, enzyme)
                cuts = cut.search(Seq(seq))
                if cuts:
                    st.success(f"{enzyme} schneidet bei {cuts}")
                else:
                    st.info("Keine Schnittstelle gefunden.")
            except Exception as e:
                st.error(f"Fehler: {e}")

    # ---------- Plasmid Map ----------
    with tabs[4]:
        st.subheader("Plasmid-Map Visualisierung")
        seq = st.text_area("DNA-Sequenz oder FASTA", height=120, key="pm_seq")
        circular = st.checkbox("Zirkulär darstellen", True)
        if st.button("Plasmid anzeigen"):
            dna = seq
            if ">" in seq:  # FASTA heuristik
                recs = list(SeqIO.parse(io.StringIO(seq), "fasta"))
                dna = str(recs[0].seq) if recs else ""
            if dna:
                length = len(dna)
                feats = [(0, length, "Backbone")]
                fig = circular_map(length, feats) if circular else linear_feature_plot(length, feats)
                st.plotly_chart(fig, use_container_width=True)
            else:
                st.info("Keine gültige DNA gefunden.")

    # ---------- Plasmid Annotator ----------
    with tabs[5]:
        st.subheader("Plasmid Annotator")
        gb_up = st.file_uploader("GenBank/FASTA laden", type=["gb","gbk","genbank","fasta","fa"], key="ann_up")
        if gb_up:
            text = gb_up.getvalue().decode()
            records = []
            if gb_up.name.lower().endswith((".gb", ".gbk", ".genbank")):
                records = list(SeqIO.parse(io.StringIO(text), "genbank"))
            else:
                records = list(SeqIO.parse(io.StringIO(text), "fasta"))
            if not records:
                st.warning("Keine Sequenz erkannt.")
            else:
                rec = records[0]
                seq = str(rec.seq).upper()
                found = []
                # vorhandene GenBank-Features mitnehmen
                if getattr(rec, "features", None):
                    for f in rec.features:
                        if isinstance(f.location, FeatureLocation):
                            found.append((int(f.location.start), int(f.location.end), f.type))
                # heuristische Motive
                for name, pats in MOTIFS.items():
                    for p in pats:
                        for s, e in search_motif(seq, p):
                            found.append((s, e, name))
                # Darstellung
                if not found:
                    st.info("Keine Motive/Features erkannt.")
                else:
                    found.sort(key=lambda x: x[0])
                    df = pd.DataFrame(found, columns=["Start","Ende","Feature"])
                    st.dataframe(df, use_container_width=True)
                    fig = circular_map(len(seq), found[:20])  # viele beschneiden
                    st.plotly_chart(fig, use_container_width=True)

    # ---------- Cloning-Primer Generator ----------
    with tabs[6]:
        st.subheader("Cloning-Primer Generator")
        seq = st.text_area("Template-Sequenz (5'→3')", height=120, key="cp_seq")
        left_overhang = st.text_input("5' Overhang Left (z. B. GAATTC)", "")
        right_overhang = st.text_input("5' Overhang Right (z. B. TCTAGA)", "")
        opt_len = st.number_input("opt. Primerlänge (bp)", 14, 35, 20)
        tm_target = st.number_input("Ziel-Tm (°C)", 50, 70, 60)
        if st.button("Primer vorschlagen"):
            if not P3_OK:
                st.error("primer3 nicht installiert. Bitte 'primer3-py' in requirements aufnehmen.")
            elif not seq:
                st.info("Template eingeben.")
            else:
                # minimaler primer3-call (wir nutzen nur left/right SINGLE primer)
                res = primer3.bindings.designPrimers(
                    {"SEQUENCE_TEMPLATE": seq},
                    {
                        "PRIMER_OPT_SIZE": int(opt_len),
                        "PRIMER_MIN_SIZE": 14,
                        "PRIMER_MAX_SIZE": 35,
                        "PRIMER_OPT_TM": float(tm_target),
                        "PRIMER_MIN_TM": float(tm_target-5),
                        "PRIMER_MAX_TM": float(tm_target+5),
                        "PRIMER_NUM_RETURN": 1
                    }
                )
                L = res.get("PRIMER_LEFT_0_SEQUENCE", "")
                R = res.get("PRIMER_RIGHT_0_SEQUENCE", "")
                if L and R:
                    L_full = (left_overhang.upper() + L)
                    R_full = (right_overhang.upper() + R)
                    st.success("Primer generiert:")
                    st.write(f"Left:  {L_full}  (GC {gc(L_full):.1f} %)")
                    st.write(f"Right: {R_full}  (GC {gc(R_full):.1f} %)")
                    fasta = f">Left\n{L_full}\n>Right\n{R_full}\n"
                    st.download_button("FASTA exportieren", fasta.encode("utf-8"),
                                       file_name="cloning_primers.fasta", mime="text/plain")
                else:
                    st.warning("Kein Primer gefunden – Parameter anpassen.")

    # ---------- Gateway / TA / TOPO ----------
    with tabs[7]:
        st.subheader("Gateway / TA / TOPO – Checks")
        seq = st.text_area("Insert-Sequenz", height=120, key="gt_seq")
        # Gateway attB-Motive (kurz)
        attb_motifs = {"attB1": "GGGGACAAGTTTGTACAAAAAAGCAGGCT",
                       "attB2": "GGGGACCACTTTGTACAAGAAAGCTGGGT"}
        if st.button("Analysieren"):
            if not seq:
                st.info("Sequenz eingeben.")
            else:
                msgs = []
                s = seq.upper()
                # TA-Klonierung: A-overhang-kompatible PCR-Produkte → oft A am 3'-Ende
                if s.endswith("A"):
                    msgs.append("TA-compatible: Insert endet mit 'A'.")
                else:
                    msgs.append("TA-Hinweis: Insert endet nicht mit 'A'.")
                # TOPO (blunt/TA) – hier einfache Hinweise
                if s.startswith("CACC"):
                    msgs.append("TOPO (Directional) Hinweis: 5'-CACC gefunden.")
                # Gateway
                for name, motif in attb_motifs.items():
                    if motif in s:
                        msgs.append(f"Gateway: {name}-Motiv erkannt.")
                st.write("\n".join(f"• {m}" for m in msgs))

    # ---------- GenBank-Viewer ----------
    with tabs[8]:
        st.subheader("GenBank-Datei Viewer")
        up = st.file_uploader("GenBank laden", type=["gb","gbk","genbank"], key="gbv_up")
        if up:
            text = up.getvalue().decode()
            recs = list(SeqIO.parse(io.StringIO(text), "genbank"))
            if not recs:
                st.warning("Keine GenBank-Features erkannt.")
            else:
                rec = recs[0]
                seq = str(rec.seq)
                feats = []
                for f in rec.features:
                    if isinstance(f.location, FeatureLocation):
                        s = int(f.location.start)
                        e = int(f.location.end)
                        feats.append((s, e, f.type))
                st.write(f"Features: {len(feats)} · Länge: {len(seq)} bp")
                df = pd.DataFrame(feats, columns=["Start","Ende","Typ"]).sort_values("Start")
                st.dataframe(df, use_container_width=True)
                fig = circular_map(len(seq), feats[:30])
                st.plotly_chart(fig, use_container_width=True)
