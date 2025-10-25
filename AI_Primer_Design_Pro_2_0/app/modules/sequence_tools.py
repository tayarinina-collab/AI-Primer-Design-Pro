from Bio.Seq import Seq
from Bio.SeqUtils import GC, MeltingTemp as mt
import re
import pandas as pd

def detect_sequence_type(seq):
    seq = seq.upper().replace("\n", "")
    if re.fullmatch(r"[ACGT]+", seq):
        return "DNA"
    elif re.fullmatch(r"[ACGU]+", seq):
        return "RNA"
    elif re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+", seq):
        return "Protein"
    else:
        return "Unknown"

def compute_basic_properties(seq):
    seq_type = detect_sequence_type(seq)
    seq_obj = Seq(seq)
    length = len(seq)
    gc_content = GC(seq) if seq_type in ["DNA", "RNA"] else None
    tm = mt.Tm_Wallace(seq) if seq_type == "DNA" else None
    return {"Type": seq_type, "Length": length, "GC%": gc_content, "Tm": tm}

def find_orfs(seq):
    seq = seq.upper()
    start_codon = "ATG"
    stop_codons = ["TAA", "TAG", "TGA"]
    orfs = []
    for frame in range(3):
        start = None
        for i in range(frame, len(seq), 3):
            codon = seq[i:i+3]
            if codon == start_codon and start is None:
                start = i
            elif codon in stop_codons and start is not None:
                orfs.append((start, i+3))
                start = None
    return orfs

def find_motifs(seq):
    motifs = {
        "TATA-box": r"TATA[AT]A[AT][AG]",
        "EcoRI": r"GAATTC",
        "BamHI": r"GGATCC",
        "HindIII": r"AAGCTT"
    }
    found = {}
    for name, pattern in motifs.items():
        found[name] = [m.start() for m in re.finditer(pattern, seq.upper())]
    return found

def gc_profile(seq, window=50):
    seq = seq.upper()
    values = []
    for i in range(0, len(seq)-window+1, window):
        window_seq = seq[i:i+window]
        gc = GC(window_seq)
        values.append(gc)
    df = pd.DataFrame({"Window": range(1, len(values)+1), "GC%": values})
    return df
