# -*- coding: utf-8 -*-
"""
Protein Tools Modul
• Translation / Reverse Translation
• Hydrophobizität (Kyte-Doolittle), pI, Molekulargewicht
• Motif-Suche (PROSITE, Pfam-ähnlich)
• Sekundärstruktur-Schätzung (rudimentär)
• 3D-Struktur-Viewer (Cloud-kompatibel)
• Protein-Domain Annotation mit UniProt (tabellarisch)
"""
import io, requests
import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Bio.Seq import Seq
from Bio.SeqUtils import ProtParam

# optional 3D viewer
try:
    import py3Dmol
    PDB_OK = True
except Exception:
    PDB_OK = False


# ---------------- Helper ----------------------------------------------------

HYDRO_SCALE = {
    'A': 1.8, 'C': 2.5, 'D': -3.5, 'E': -3.5, 'F': 2.8, 'G': -0.4, 'H': -3.2, 'I': 4.5,
    'K': -3.9, 'L': 3.8, 'M': 1.9, 'N': -3.5, 'P': -1.6, 'Q': -3.5, 'R': -4.5,
    'S': -0.8, 'T': -0.7, 'V': 4.2, 'W': -0.9, 'Y': -1.3
}

def hydrophobicity_profile(seq: str, window: int = 9):
    vals = [HYDRO_SCALE.get(a, 0) for a in seq]
    smoothed = []
    for i in range(len(vals)):
        w = vals[max(0, i - window // 2):min(len(vals), i + window // 2)]
        smoothed.append(np.mean(w))
    return smoothed

def simple_secondary_structure(seq: str):
    """rudimentäre Klassifikation auf Basis von Aminosäurenhäufigkeit"""
    alpha_pref = "MALEK"
    beta_pref = "VIYFW"
    coil_pref = "PGNSTQ"
    if not seq:
        return {}
    comp = {k: seq.count(k) / len(seq) for k in "ACDEFGHIKLMNPQRSTVWY"}
    alpha = sum(comp.get(a, 0) for a in alpha_pref)
    beta = sum(comp.get(a, 0) for a in beta_pref)
    coil = sum(comp.get(a, 0) for a in coil_pref)
    total = alpha + beta + coil
    if total == 0:
        total = 1
    return {"α-Helix": alpha / total, "β-Strang": beta / total, "coil": coil / total}

def reverse_translate(protein: str):
    """reverse translation (codon usage simplified for E. coli)"""
    codons = {
        'A': "GCT", 'C': "TGC", 'D': "GAT", 'E': "GAA", 'F': "TTT", 'G': "GGT", 'H': "CAT",
        'I': "ATT", 'K': "AAA", 'L': "CTG", 'M': "ATG", 'N': "AAT", 'P': "CCT", 'Q': "CAA",
        'R': "CGT", 'S': "TCT", 'T': "ACC", 'V': "GTG", 'W': "TGG", 'Y': "TAT"
    }
    return "".join(codons.get(a, "NNN") for a in protein)

def motif_search(seq: str):
    """rudimentäre PROSITE-/Pfam-Simulation: erkennt häufige Motive"""
    motifs = {
        "N-glycosylation": r"N[^P][ST][^P]",
        "Protein kinase ATP binding": r"[LIV]-x(2)-[LIV]-x(3)-G-[E]",
        "Zinc finger C2H2": r"C.{2,4}C.{12}H.{3,5}H",
        "Leucine zipper": r"L.{6}L.{6}L.{6}L",
        "Nuclear localization signal": r"K.{2,3}K",
    }
    results = []
    import re
    for name, pattern in motifs.items():
        hits = [(m.start(), m.end()) for m in re.finditer(pattern, seq)]
        if hits:
            results.append((name, len(hits), hits))
    return results


# ---------------- UI -------------------------------------------------------

def run_protein_tools():
    st.title("🧬 Protein Tools")
    st.caption("Translation · Hydrophobizität · Motif-Suche · Struktur & 3D-Viewer · UniProt Annotation")

    tabs = st.tabs([
        "Translation / Reverse",
        "Hydrophobizität & pI",
        "Motif-Suche",
        "Sekundärstruktur",
        "3D-Struktur",
        "UniProt Annotation"
    ])

    # ---------- Translation ----------
    with tabs[0]:
        st.subheader("DNA↔Protein Übersetzung")
        dna = st.text_area("DNA-Sequenz (5'→3')", height=120, key="dna_seq")
        if st.button("Übersetzen"):
            seq = Seq(dna.replace("U", "T"))
            prot = seq.translate(to_stop=True)
            st.success(f"Protein ({len(prot)} aa):")
            st.code(str(prot))
        st.markdown("---")
        protein = st.text_area("Proteinsequenz", height=120, key="prot_seq")
        if st.button("Reverse Translation"):
            rt = reverse_translate(protein)
            st.success("Rückübersetzung (E. coli Codon Usage):")
            st.code(rt)

    # ---------- Hydrophobizität ----------
    with tabs[1]:
        st.subheader("Hydrophobizität / Molekulargewicht / pI")
        prot = st.text_area("Proteinsequenz (1-letter code)", height=120, key="hydro_seq")
        if st.button("Analysieren", key="analyze_hydro"):
            if prot:
                X = ProtParam.ProteinAnalysis(prot)
                mw = X.molecular_weight()
                pi = X.isoelectric_point()
                hyd = hydrophobicity_profile(prot)
                st.write(f"**Molekulargewicht:** {mw / 1000:.2f} kDa  ·  **pI:** {pi:.2f}")
                fig, ax = plt.subplots()
                ax.plot(hyd)
                ax.set_title("Kyte-Doolittle Hydrophobizität")
                ax.set_xlabel("Position")
                ax.set_ylabel("Hydrophobizität")
                st.pyplot(fig, use_container_width=True)

    # ---------- Motif-Suche ----------
    with tabs[2]:
        st.subheader("Motif-Suche (PROSITE/Pfam-ähnlich)")
        prot = st.text_area("Proteinsequenz", height=120, key="motif_seq")
        if st.button("Motive finden"):
            results = motif_search(prot)
            if not results:
                st.info("Keine Standardmotive erkannt.")
            else:
                df = pd.DataFrame(
                    [{"Motif": n, "Treffer": c, "Positionen": str(p)} for n, c, p in results]
                )
                st.dataframe(df, use_container_width=True)

    # ---------- Sekundärstruktur ----------
    with tabs[3]:
        st.subheader("Sekundärstruktur-Vorhersage (rudimentär)")
        prot = st.text_area("Proteinsequenz", height=120, key="ss_seq")
        if st.button("Vorhersage starten"):
            pred = simple_secondary_structure(prot)
            if not pred:
                st.warning("Bitte eine gültige Sequenz eingeben.")
            else:
                st.write(pred)
                fig, ax = plt.subplots()
                ax.bar(pred.keys(), pred.values(), color=["#64b5f6", "#81c784", "#ffb74d"])
                ax.set_ylabel("Anteil")
                st.pyplot(fig, use_container_width=True)

    # ---------- 3D-Struktur ----------
    with tabs[4]:
        st.subheader("3D-Struktur-Viewer (PDB) 🌐")
        pdb_id = st.text_input("PDB ID (z. B. 1CRN oder 6M0J)", "")
        if st.button("Struktur anzeigen"):
            import streamlit.components.v1 as components
            if not pdb_id:
                st.info("Bitte PDB-ID eingeben.")
            else:
                try:
                    pdb_embed_url = f"https://www.rcsb.org/3d-view/{pdb_id}"
                    st.components.v1.iframe(pdb_embed_url, height=600, width=800)
                except Exception as e:
                    st.error(f"Fehler beim Laden der PDB-Ansicht: {e}")

    # ---------- UniProt Annotation ----------
    with tabs[5]:
        st.subheader("Protein-Domain Annotation (UniProt)")
        uniprot_id = st.text_input("UniProt ID (z. B. P69905 für HBA_HUMAN)", "")
        if st.button("Daten abrufen"):
            if not uniprot_id:
                st.info("Bitte ID eingeben.")
            else:
                url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.txt"
                try:
                    r = requests.get(url, timeout=10)
                    if r.status_code == 200:
                        raw_text = r.text
                        features = [l for l in raw_text.splitlines() if l.startswith("FT")]

                        def parse_uniprot_features(lines):
                            rows = []
                            for line in lines:
                                parts = line.split()
                                if len(parts) > 2:
                                    feature = parts[1]
                                    rest = " ".join(parts[2:])
                                    rows.append({"Feature": feature, "Beschreibung": rest})
                            return pd.DataFrame(rows)

                        df = parse_uniprot_features(features)

                        st.success(f"Eintrag erfolgreich geladen: **{uniprot_id}**  ·  {len(df)} Features erkannt")
                        st.dataframe(df, use_container_width=True)

                        with st.expander("Rohdaten anzeigen (Original UniProt Format)"):
                            st.text_area("UniProt Raw Feature Data", "\n".join(features), height=300)
                    else:
                        st.warning("Kein Eintrag gefunden oder ungültige ID.")
                except Exception as e:
                    st.error(f"Fehler beim Abrufen: {e}")
