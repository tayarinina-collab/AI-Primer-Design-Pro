# -*- coding: utf-8 -*-
"""
Plasmid Designer – Visualisierung, Import, Restriktionskarte, Kloning & Export
Offline-fähig; nutzt Biopython (SeqIO), matplotlib, pandas.
"""

from __future__ import annotations
import io
import math
from dataclasses import dataclass
from typing import List, Tuple, Dict

import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

# --- Optional: Bio.Restriction, wenn vorhanden (sonst Fallback) -----
try:
    from Bio.Restriction import RestrictionBatch, AllEnzymes
    BIO_RESTRICTION_OK = True
except Exception:
    BIO_RESTRICTION_OK = False


# =========================== Datenstrukturen =================================

@dataclass
class Feature:
    name: str
    start: int
    end: int
    ftype: str = "misc_feature"
    strand: int = 1


FEATURE_COLORS = {
    "ori":        "#FFAC33",
    "promoter":   "#70C1B3",
    "cds":        "#C16E70",
    "terminator": "#F25F5C",
    "marker":     "#247BA0",
    "misc":       "#9CA3AF",
}


# =========================== Hilfsfunktionen =================================

def wrap_len(idx: int, L: int) -> int:
    return idx % L


def normalize_feature(start: int, end: int, L: int) -> List[Tuple[int, int]]:
    start, end = wrap_len(start, L), wrap_len(end, L)
    if start < end:
        return [(start, end)]
    elif start > end:
        return [(start, L), (0, end)]
    else:
        return []


def gc_percent(seq: str) -> float:
    s = seq.upper().replace("U", "T")
    return 100.0 * (s.count("G") + s.count("C")) / max(1, len(s))


def revcomp(seq: str) -> str:
    return str(Seq(seq).reverse_complement())


def find_orfs(seq: str, min_aa: int = 100) -> List[Tuple[int, int]]:
    s = seq.upper()
    stops = {"TAA", "TAG", "TGA"}
    L = len(s)
    orfs = []
    for frame in range(3):
        i = frame
        while i + 3 <= L:
            cod = s[i:i+3]
            if cod == "ATG":
                j = i + 3
                while j + 3 <= L:
                    if s[j:j+3] in stops:
                        aa_len = (j - i) // 3
                        if aa_len >= min_aa:
                            orfs.append((i, j+3))
                        i = j
                        break
                    j += 3
            i += 3
    return orfs


BASIC_ENZYMES = {
    "EcoRI": "GAATTC",
    "BamHI": "GGATCC",
    "XbaI": "TCTAGA",
    "HindIII": "AAGCTT",
    "PstI": "CTGCAG",
    "KpnI": "GGTACC",
}


def find_sites_simple(seq: str, motif: str) -> List[int]:
    s = seq.upper()
    m = motif.upper()
    rc = revcomp(m)
    L = len(s)
    hits = []
    k = len(m)
    for i in range(L - k + 1):
        if s[i:i + k] == m or s[i:i + k] == rc:
            hits.append(i)
    return hits


def make_seqrecord(seq: str, name="plasmid", description="") -> SeqRecord:
    return SeqRecord(Seq(seq), id=name, name=name, description=description)


# =========================== Darstellung =====================================

def draw_plasmid(seq: str, features: List[Feature], title: str = "Plasmid-Karte"):
    L = len(seq)
    if L == 0:
        st.warning("Keine Sequenz vorhanden.")
        return

    fig, ax = plt.subplots(figsize=(7, 7), subplot_kw={'projection': 'polar'})
    ax.set_title(title)
    ax.set_theta_direction(-1)
    ax.set_theta_offset(math.pi / 2.0)

    # Grundkreis
    ax.plot(np.linspace(0, 2 * math.pi, 360), [1.0] * 360, lw=2, color="#555")

    # Tick Marks (alle 1000 bp)
    ticks = max(1, L // 1000)
    for t in range(ticks + 1):
        pos = int(t * L / ticks)
        theta = (pos / L) * 2 * math.pi
        ax.plot([theta, theta], [0.95, 1.05], color="#999", lw=1)
        ax.text(theta, 1.1, f"{pos}", fontsize=8,
                rotation=-(theta * 180 / math.pi - 90), ha="center", va="center")

    # Features zeichnen
    for f in features:
        color = FEATURE_COLORS.get(f.ftype, FEATURE_COLORS["misc"])
        parts = normalize_feature(f.start, f.end, L)
        for (a, b) in parts:
            th0 = (a / L) * 2 * math.pi
            th1 = (b / L) * 2 * math.pi
            ax.bar(x=(th0 + th1) / 2, height=0.08, width=th1 - th0,
                   bottom=0.9, color=color, edgecolor="#333", alpha=0.7)
            mid = (a + (b - a) / 2) % L
            tm = (mid / L) * 2 * math.pi
            ax.text(tm, 0.85, f.name, fontsize=8,
                    rotation=-(tm * 180 / math.pi - 90), ha="center", va="center")

    ax.set_yticklabels([])
    ax.set_xticklabels([])
    ax.set_ylim(0.7, 1.2)
    st.pyplot(fig, use_container_width=True)
    return fig


def fig_to_svg_bytes(fig) -> bytes:
    buf = io.BytesIO()
    fig.savefig(buf, format="svg", bbox_inches="tight")
    return buf.getvalue()


# =========================== Haupt-UI ========================================

def run_plasmid_designer():
    st.title("🧫 Plasmid-Karte (Designer)")

    if "plasmid_seq" not in st.session_state:
        st.session_state.plasmid_seq = ""
    if "plasmid_features" not in st.session_state:
        st.session_state.plasmid_features: List[Feature] = []

    tabs = st.tabs([
        "🧬 Plasmid-Karte",
        "📄 GenBank/FASTA Import",
        "✂️ Restriktionskarte",
        "🧪 Kloning-Simulation",
        "📤 Export",
    ])

    # -------------------------------------------------------------------------
    # TAB 1 – Plasmid-Karte
    # -------------------------------------------------------------------------
    with tabs[0]:
        st.subheader("Plasmid-Karte & Feature-Management")

        seq = st.text_area("Sequenz (DNA, 5'→3')",
                           value=st.session_state.plasmid_seq,
                           height=160)

        if seq.strip():
            length = len(seq)
            gc = gc_percent(seq)
            st.caption(f"Länge: **{length} bp**, GC: **{gc:.1f}%**")

            if not st.session_state.plasmid_features:
                st.session_state.plasmid_features = [
                    Feature("ori", 100, 300, "ori", 1),
                    Feature("AmpR", 350, 550, "marker", 1),
                    Feature("MCS", 560, 600, "misc", 1),
                ]

            df = pd.DataFrame([f.__dict__ for f in st.session_state.plasmid_features])
            edited = st.data_editor(df, num_rows="dynamic", use_container_width=True)

            if st.button("Änderungen übernehmen / Update"):
                feats = []
                for _, r in edited.iterrows():
                    try:
                        feats.append(Feature(str(r["name"]), int(r["start"]), int(r["end"]),
                                             str(r["ftype"]), int(r["strand"])))
                    except Exception:
                        pass
                st.session_state.plasmid_features = feats
                st.session_state.plasmid_seq = seq
                st.success("Aktualisiert ✅")

            draw_plasmid(seq, st.session_state.plasmid_features)

        else:
            st.info("🔹 Bitte zuerst eine Sequenz eingeben oder im Tab **GenBank/FASTA Import** laden.")

    # -------------------------------------------------------------------------
    # TAB 2 – GenBank / FASTA Import
    # -------------------------------------------------------------------------
    with tabs[1]:
        st.subheader("GenBank / FASTA Import")
        upl = st.file_uploader("Datei wählen", type=["gb", "gbk", "genbank", "fasta", "fa", "txt"])
        if upl is not None:
            raw = upl.read().decode("utf-8", errors="ignore")
            fmt = "genbank" if upl.name.lower().endswith(("gb", "gbk", "genbank")) else "fasta"
            try:
                recs = list(SeqIO.parse(io.StringIO(raw), fmt))
                if not recs:
                    st.warning("Keine Einträge gefunden.")
                else:
                    rec = recs[0]
                    st.session_state.plasmid_seq = str(rec.seq)
                    feats = []
                    if fmt == "genbank" and hasattr(rec, "features"):
                        for f in rec.features:
                            try:
                                start = int(f.location.start)
                                end = int(f.location.end)
                                strand = int(f.location.strand or 1)
                                name = f.qualifiers.get("label", ["feature"])[0]
                                ftype = f.type or "misc"
                                feats.append(Feature(name, start, end, ftype, strand))
                            except Exception:
                                pass
                    st.session_state.plasmid_features = feats
                    st.success(f"Importiert ✅ | Länge: {len(rec.seq)} bp | Features: {len(feats)}")
            except Exception as e:
                st.error(f"Fehler beim Import: {e}")

    # -------------------------------------------------------------------------
    # TAB 3 – Restriktionskarte
    # -------------------------------------------------------------------------
    with tabs[2]:
        st.subheader("✂️ Restriktionskarte")

        seq = st.session_state.get("plasmid_seq", "").upper().replace("\n", "").replace(" ", "")
        if not seq:
            st.info("⚠️ Keine Sequenz gefunden. Bitte zuerst eine DNA-Sequenz eingeben oder importieren.")
        else:
            st.caption(f"Länge: **{len(seq)} bp**, GC: **{gc_percent(seq):.1f}%**")

            enzyme_list = list(BASIC_ENZYMES.keys())
            picked = st.multiselect("Restriktionsenzyme auswählen",
                                    options=enzyme_list,
                                    default=["EcoRI", "BamHI", "HindIII"])

            if st.button("🔎 Restriktionsschnittstellen suchen"):
                rows = []
                for name in picked:
                    motif = BASIC_ENZYMES[name]
                    hits = find_sites_simple(seq, motif)
                    for pos in hits:
                        rows.append({"Enzym": name, "Position (bp)": pos})
                if rows:
                    df_sites = pd.DataFrame(rows).sort_values("Position (bp)")
                    st.success(f"{len(rows)} Schnittstellen gefunden ✅")
                    st.dataframe(df_sites, use_container_width=True)

                    fig, ax = plt.subplots(figsize=(6, 1))
                    ax.set_xlim(0, len(seq))
                    ax.set_ylim(0, 1)
                    ax.set_yticks([])
                    for _, r in df_sites.iterrows():
                        ax.axvline(x=r["Position (bp)"], color="red", lw=1)
                        ax.text(r["Position (bp)"], 0.5, r["Enzym"], rotation=90,
                                va="center", ha="center", fontsize=8, color="red")
                    ax.set_xlabel("Position (bp)")
                    ax.set_title("Restriktionsschnittstellen")
                    st.pyplot(fig, use_container_width
