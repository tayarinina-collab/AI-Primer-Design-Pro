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

# --- Optional: Bio.Restriction, wenn vorhanden (sonst eigene Erkennung) -----
try:
    from Bio.Restriction import RestrictionBatch, AllEnzymes
    BIO_RESTRICTION_OK = True
except Exception:
    BIO_RESTRICTION_OK = False


# =========================== Datenstrukturen =================================

@dataclass
class Feature:
    name: str
    start: int  # 0-based, inclusive
    end: int    # 0-based, exclusive
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
    """Begrenze Index in [0, L)."""
    idx %= L
    return idx


def normalize_feature(start: int, end: int, L: int) -> List[Tuple[int, int]]:
    """
    Normalisiere ein Feature (start,end) auf 0..L-1. Liefert eine oder zwei
    nicht-überlappende Intervalle, falls das Feature über den Origin läuft.
    """
    start, end = wrap_len(start, L), wrap_len(end, L)
    if start < end:
        return [(start, end)]
    elif start > end:
        return [(start, L), (0, end)]
    else:
        # volle Runde oder null – hier Null-Länge ignorieren
        return []


def gc_percent(seq: str) -> float:
    s = seq.upper().replace("U", "T")
    return 100.0 * (s.count("G") + s.count("C")) / max(1, len(s))


def revcomp(seq: str) -> str:
    return str(Seq(seq).reverse_complement())


def find_orfs(seq: str, min_aa: int = 100) -> List[Tuple[int, int]]:
    """Einfache ORF-Erkennung (ATG...*), returns ORFs als (start,end) in nt."""
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


# Fallback Restriktionsenzyme (falls Bio.Restriction nicht verfügbar ist)
BASIC_ENZYMES = {
    "EcoRI":  "GAATTC",
    "BamHI":  "GGATCC",
    "XbaI":   "TCTAGA",
    "HindIII":"AAGCTT",
    "PstI":   "CTGCAG",
    "KpnI":   "GGTACC",
}


def find_sites_simple(seq: str, motif: str) -> List[int]:
    """Finde alle Startpositionen eines Motifs (und Reverse-Complement) in seq."""
    s = seq.upper()
    m = motif.upper()
    rc = revcomp(m)
    L = len(s)
    hits = []
    k = len(m)
    for i in range(L - k + 1):
        if s[i:i+k] == m or s[i:i+k] == rc:
            hits.append(i)
    return hits


def make_seqrecord(seq: str, name="plasmid", description="") -> SeqRecord:
    rec = SeqRecord(Seq(seq), id=name, name=name, description=description)
    return rec


# =========================== Darstellung =====================================

def draw_plasmid(seq: str, features: List[Feature], title: str = "Plasmid-Karte"):
    """
    Kreisplot der Features. Nutzt einfache Polarkoordinaten.
    """
    L = len(seq)
    if L == 0:
        st.warning("Keine Sequenz vorhanden.")
        return

    fig, ax = plt.subplots(figsize=(7, 7), subplot_kw={'projection': 'polar'})
    ax.set_title(title)
    ax.set_theta_direction(-1)  # clockwise
    ax.set_theta_offset(math.pi / 2.0)  # start at top (12 o'clock)

    # Grundkreis
    ax.plot(np.linspace(0, 2*math.pi, 360), [1.0]*360, lw=2, color="#555")

    # Tick Marks (alle 1000 bp)
    ticks = max(1, L // 1000)
    for t in range(ticks+1):
        pos = int(t * L / ticks)
        theta = (pos / L) * 2 * math.pi
        ax.plot([theta, theta], [0.95, 1.05], color="#999", lw=1)
        ax.text(theta, 1.1, f"{pos}", fontsize=8,
                rotation=-(theta * 180 / math.pi - 90), ha="center", va="center")

    # Features
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

    # --- Session State für Sequenz & Features ---
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

    # =========================================================================
    # TAB 1 – Plasmid-Karte
    # =========================================================================
    with tabs[0]:
        st.subheader("Plasmid-Karte & Feature-Management")

        # Sequenz Eingabe
        seq = st.text_area(
            "Sequenz (DNA, 5'→3')",
            value=st.session_state.plasmid_seq,
            height=140,
            help="Hier die Vektor- oder Plasmidsequenz einfügen.",
        )

        # Info: Länge & GC
        if seq.strip():
            st.caption(f"Länge: **{len(seq)} bp**, GC: **{gc_percent(seq):.1f}%**")

        # Feature-Tabelle
        st.markdown("**Features (bearbeitbar)**")
        if not st.session_state.plasmid_features and seq.strip():
            # kleinem Beispiel, wenn frisch
            st.session_state.plasmid_features = [
                Feature("ori", 100, 300, "ori", 1),
                Feature("AmpR", 800, 1600, "marker", 1),
                Feature("MCS", 1600, 1800, "misc", 1),
            ]

        # Editierbare Tabelle
        df = pd.DataFrame(
            [
                {
                    "name": f.name,
                    "start": f.start,
                    "end": f.end,
                    "type": f.ftype,
                    "strand": f.strand,
                }
                for f in st.session_state.plasmid_features
            ]
        )
        edited = st.data_editor(
            df,
            num_rows="dynamic",
            use_container_width=True,
            key="features_editor",
            column_config={
                "type": st.column_config.SelectboxColumn(
                    "type",
                    options=list(FEATURE_COLORS.keys()) + ["misc"],
                ),
                "strand": st.column_config.NumberColumn(
                    "strand", min_value=-1, max_value=1, step=1
                ),
            },
        )

        # Übernehmen-Button
        if st.button("Änderungen übernehmen / Update"):
            feats: List[Feature] = []
            for _, r in edited.iterrows():
                try:
                    feats.append(
                        Feature(
                            name=str(r["name"]),
                            start=int(r["start"]),
                            end=int(r["end"]),
                            ftype=str(r["type"]) if pd.notna(r["type"]) else "misc",
                            strand=int(r["strand"]) if pd.notna(r["strand"]) else 1,
                        )
                    )
                except Exception:
                    pass
            st.session_state.plasmid_features = feats
            st.session_state.plasmid_seq = seq
            st.success("Aktualisiert ✅")

        # Plot
        if seq.strip():
            fig = draw_plasmid(seq, st.session_state.plasmid_features, "Plasmid-Karte")
        else:
            st.info("🔹 Bitte zuerst eine Sequenz eingeben oder in Tab **GenBank/FASTA Import** laden.")

    # =========================================================================
    # TAB 2 – GenBank / FASTA Import
    # =========================================================================
    with tabs[1]:
        st.subheader("GenBank / FASTA Import")

        upl = st.file_uploader("Datei wählen", type=["gb", "gbk", "genbank", "fasta", "fa", "txt"])
        if upl is not None:
            raw = upl.read().decode("utf-8", errors="ignore")
            fmt = "genbank" if upl.name.lower().endswith(("gb","gbk","genbank")) else "fasta"
            try:
                recs = list(SeqIO.parse(io.StringIO(raw), fmt))
                if not recs:
                    st.warning("Keine Einträge gefunden.")
                else:
                    rec = recs[0]
                    st.session_state.plasmid_seq = str(rec.seq)
                    feats: List[Feature] = []
                    if fmt == "genbank" and hasattr(rec, "features"):
                        for f in rec.features:
                            try:
                                start = int(f.location.nofuzzy_start)
                                end = int(f.location.nofuzzy_end)
                                strand = int(f.location.strand or 1)
                                name = f.qualifiers.get("label", f.qualifiers.get("gene", f.qualifiers.get("note", ["feature"])))[0]
                                ftype = str(f.type or "misc")
                                feats.append(Feature(name, start, end, ftype, strand))
                            except Exception:
                                pass
                    else:
                        # fallback: einfache ORFs vorschlagen
                        for (a,b) in find_orfs(str(rec.seq), 100):
                            feats.append(Feature(f"ORF_{a}", a, b, "cds", 1))

                    st.session_state.plasmid_features = feats
                    st.success(f"Importiert ✅  | Länge: {len(rec.seq)} bp | Features: {len(feats)}")
            except Exception as e:
                st.error(f"Fehler beim Import: {e}")

    # =========================================================================
    # TAB 3 – Restriktionskarte
    # =========================================================================
    with tabs[2]:
        st.subheader("Restriktionskarte")
        if not st.session_state.plasmid_seq:
            st.info("Bitte zuerst eine Sequenz eingeben (Tab **Plasmid-Karte**) oder importieren.")
        else:
            seq = st.session_state.plasmid_seq
            L = len(seq)

            if BIO_RESTRICTION_OK:
                st.caption("🔧 Analyse mit **Bio.Restriction** (falls ein Enzym-Set ausgewählt wird, kann es etwas dauern)…")
                enzyme_names = st.multiselect(
                    "Enzyme auswählen (optional – leer = kleines Standard-Set)",
                    options=[e.__name__ for e in AllEnzymes],
                    default=[]
                )
                if st.button("Schnittstellen berechnen"):
                    if enzyme_names:
                        rb = RestrictionBatch(enzyme_names)
                    else:
                        rb = RestrictionBatch(list(BASIC_ENZYMES.keys()))
                    analysis = rb.search(Seq(seq))
                    rows = []
                    for enz, hits in analysis.items():
                        for h in hits:
                            rows.append({"Enzym": enz, "Position": int(h)})
                    df_sites = pd.DataFrame(rows).sort_values("Position")
                    if df_sites.empty:
                        st.info("Keine Schnittstellen gefunden.")
                    else:
                        st.dataframe(df_sites, use_container_width=True)
            else:
                st.caption("🧰 Einfacher Fallback-Scanner (ohne Bio.Restriction).")
                picked = st.multiselect(
                    "Standard-Enzyme",
                    options=list(BASIC_ENZYMES.keys()),
                    default=["EcoRI", "BamHI", "XbaI"]
                )
                if st.button("Schnittstellen scannen (Fallback)"):
                    rows = []
                    for name in picked:
                        motif = BASIC_ENZYMES[name]
                        for pos in find_sites_simple(seq, motif):
                            rows.append({"Enzym": name, "Position": pos})
                    df_sites = pd.DataFrame(rows).sort_values("Position")
                    if df_sites.empty:
                        st.info("Keine Schnittstellen gefunden.")
                    else:
                        st.dataframe(df_sites, use_container_width=True)

    # =========================================================================
    # TAB 4 – Kloning-Simulation
    # =========================================================================
    with tabs[3]:
        st.subheader("Kloning-Simulation (einfach)")
        if not st.session_state.plasmid_seq:
            st.info("Bitte zuerst eine Sequenz eingeben/importieren.")
        else:
            seq = st.session_state.plasmid_seq
            insert = st.text_area("Insert-Sequenz", height=120, placeholder="Hier Insert eingeben…")
            pos = st.number_input("Einfügeposition (bp, 0-basiert)", min_value=0, max_value=len(seq), value=0)
            rc = st.checkbox("Insert als reverse-complement einfügen", value=False)
            if st.button("Virtuell ligieren"):
                ins = revcomp(insert) if rc else insert
                new_seq = seq[:pos] + ins + seq[pos:]
                st.success(f"Neue Länge: {len(new_seq)} bp")
                st.code(new_seq[:300] + ("..." if len(new_seq) > 300 else ""), language="text")

                # Update Button
                if st.button("Als aktuelle Plasmid-Sequenz übernehmen"):
                    st.session_state.plasmid_seq = new_seq
                    st.info("Sequenz übernommen.")

    # =========================================================================
    # TAB 5 – Export
    # =========================================================================
    with tabs[4]:
        st.subheader("Export")

        if not st.session_state.plasmid_seq:
            st.info("Keine Sequenz vorhanden.")
        else:
            seq = st.session_state.plasmid_seq
            feats = st.session_state.plasmid_features

            # Export: FASTA
            rec = make_seqrecord(seq, name="plasmid", description="exported from AI Primer Design Pro")
            fasta_bytes = io.BytesIO()
            SeqIO.write(rec, fasta_bytes, "fasta")
            st.download_button(
                "⬇️ FASTA exportieren",
                data=fasta_bytes.getvalue(),
                file_name="plasmid.fasta",
                mime="text/fasta",
            )

            # Export: Feature-CSV
            feat_df = pd.DataFrame(
                [{"name": f.name, "start": f.start, "end": f.end, "type": f.ftype, "strand": f.strand} for f in feats]
            )
            st.download_button(
                "⬇️ Features als CSV",
                data=feat_df.to_csv(index=False).encode("utf-8"),
                file_name="plasmid_features.csv",
                mime="text/csv",
            )

            # Export: SVG der Plasmid-Karte
            st.markdown("Karte als **SVG** exportieren:")
            fig = draw_plasmid(seq, feats, "Plasmid-Karte (Export-Vorschau)")
            if fig:
                svg_bytes = fig_to_svg_bytes(fig)
                st.download_button(
                    "⬇️ SVG downloaden",
                    data=svg_bytes,
                    file_name="plasmid_map.svg",
                    mime="image/svg+xml",
                )

# Ende Modul
