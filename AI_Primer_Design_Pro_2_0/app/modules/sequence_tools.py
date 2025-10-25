from Bio.SeqUtils import gc_fraction
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
import re
import pandas as pd
import matplotlib.pyplot as plt
import io

# --- 1️⃣ Basis-Berechnungen ---
def compute_basic_properties(sequence: str):
    seq = Seq(sequence.upper().replace("\n", "").replace(" ", ""))
    length = len(seq)
    gc_content = gc_fraction(seq) * 100
    tm = mt.Tm_Wallace(seq)
    return {"Länge (bp)": length, "GC (%)": round(gc_content, 2), "Tm (°C)": round(tm, 2)}

# --- 2️⃣ ORF-Erkennung ---
def find_orfs(sequence: str, min_length: int = 100):
    seq = Seq(sequence.upper())
    start_codon = "ATG"
    stop_codons = ["TAA", "TAG", "TGA"]
    orfs = []
    for frame in range(3):
        trans = seq[frame:]
        for i in range(0, len(trans) - 3, 3):
            codon = str(trans[i:i + 3])
            if codon == start_codon:
                for j in range(i + 3, len(trans) - 3, 3):
                    stop = str(trans[j:j + 3])
                    if stop in stop_codons:
                        orf_length = j + 3 - i
                        if orf_length >= min_length:
                            orfs.append((frame + 1, i, j + 3, orf_length))
                        break
    df = pd.DataFrame(orfs, columns=["Frame", "Start", "Stop", "Länge"])
    return df if not df.empty else pd.DataFrame(columns=["Frame", "Start", "Stop", "Länge"])

# --- 3️⃣ Motiv-Erkennung (z. B. TATA-Box, Restriktionsstellen) ---
def find_motifs(sequence: str):
    motifs = {
        "TATA-Box": "TATA[AT]A[AT]",
        "EcoRI": "GAATTC",
        "BamHI": "GGATCC",
        "HindIII": "AAGCTT"
    }
    results = []
    for name, pattern in motifs.items():
        for match in re.finditer(pattern, sequence.upper()):
            results.append({"Motiv": name, "Position": match.start() + 1})
    return pd.DataFrame(results) if results else pd.DataFrame(columns=["Motiv", "Position"])

# --- 4️⃣ GC-Profil ---
def gc_profile(sequence: str, window_size: int = 50):
    sequence = sequence.upper().replace("\n", "").replace(" ", "")
    gc_values = []
    for i in range(0, len(sequence) - window_size + 1, window_size):
        window = sequence[i:i + window_size]
        gc = gc_fraction(window) * 100
        gc_values.append(gc)
    return gc_values

# --- 5️⃣ GC-Plot erstellen ---
def plot_gc_profile(sequence: str):
    gc_values = gc_profile(sequence)
    fig, ax = plt.subplots(figsize=(8, 3))
    ax.plot(gc_values, color="teal")
    ax.set_title("GC-Profil")
    ax.set_xlabel("Fenster (à 50 bp)")
    ax.set_ylabel("GC (%)")
    plt.tight_layout()

    buf = io.BytesIO()
    plt.savefig(buf, format="png")
    buf.seek(0)
    return buf
