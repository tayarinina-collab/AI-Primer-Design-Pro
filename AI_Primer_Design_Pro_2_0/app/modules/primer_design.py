import streamlit as st
import primer3
from Bio.Seq import Seq
import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt
from io import StringIO

# ------------------------------------------------------
#  🔬  Primer Design & PCR Tools – Vollmodul
# ------------------------------------------------------

def run_primer_design():
    st.title("🧬 Primer Design & PCR Tools")
    st.caption("Vollintegration mit Primer3 • ΔG-Analyse • qPCR-Probe-Design • CSV/FASTA-Export")

    # --- Eingabeparameter in der Sidebar ---
    with st.sidebar:
        st.header("⚙️ Parameter einstellen")
        seq_input = st.text_area("DNA-Template (FASTA oder Text)", height=150)
        product_size = st.text_input("Produktgröße (z. B. 80-200)", "80-200")
        primer_len = st.slider("Primerlänge (optimal)", 18, 30, 20)
        tm_opt = st.slider("Tm optimal (°C)", 50, 75, 60)
        gc_min = st.number_input("GC% min", 20, 80, 40)
        gc_max = st.number_input("GC% max", 20, 80, 60)
        salt_mM = st.number_input("Na⁺-Konzentration (mM)", 1.0, 200.0, 50.0)
        mg_mM = st.number_input("Mg²⁺-Konzentration (mM)", 0.0, 5.0, 1.5)
        pick_probe = st.checkbox("qPCR-Probe (TaqMan / MB) hinzufügen", value=True)
        reporter = st.selectbox("Reporter-Fluorophor", ["FAM","HEX","ROX","Cy5","VIC"])
        quencher = st.selectbox("Quencher", ["BHQ1","BHQ2","TAMRA"])
        run_btn = st.button("🚀 Primer entwerfen")

    # --- Hilfsfunktionen ---
    def clean_seq(s):
        s = s.strip()
        if s.startswith(">"):
            s = "\n".join([l for l in s.splitlines() if not l.startswith(">")])
        return re.sub(r"[^ACGTURYKMSWBDHVNacgturykmswbdhvn]", "", s).upper()

    def dimer_hairpin_scores(fwd, rev):
        sd_f = primer3.calcHomodimer(fwd).dg
        hp_f = primer3.calcHairpin(fwd).dg
        sd_r = primer3.calcHomodimer(rev).dg
        hp_r = primer3.calcHairpin(rev).dg
        cd = primer3.calcHeterodimer(fwd, rev).dg
        return sd_f, hp_f, sd_r, hp_r, cd

    def score_pair(fwd, rev):
        gc_f = (fwd.count("G")+fwd.count("C"))/len(fwd)*100
        gc_r = (rev.count("G")+rev.count("C"))/len(rev)*100
        tm_f = primer3.calcTm(fwd)
        tm_r = primer3.calcTm(rev)
        tm_diff = abs(tm_f - tm_r)
        sd_f, hp_f, sd_r, hp_r, cd = dimer_hairpin_scores(fwd, rev)
        base = 100 - (abs(gc_f-50)/2) - (abs(gc_r-50)/2) - tm_diff
        base -= max(0, (-min(hp_f, hp_r)-3))*3
        base -= max(0, (-cd-6))*2
        return int(np.clip(base, 0, 100)), {
            "gc_f": round(gc_f,1),"gc_r":round(gc_r,1),
            "tm_f":round(tm_f,1),"tm_r":round(tm_r,1),
            "tm_diff":round(tm_diff,1),"dG_cross":round(cd,1)
        }

    def plot_dg(metrics):
        fig, ax = plt.subplots()
        bars = ["Tm_F","Tm_R","ΔTm","ΔG_cross"]
        vals = [metrics["tm_f"], metrics["tm_r"], metrics["tm_diff"], metrics["dG_cross"]]
        ax.bar(bars, vals)
        ax.set_ylabel("Wert / ΔG (kcal/mol)")
        ax.set_title("Thermodynamische Bewertung")
        return fig

    # --- Hauptfunktion ---
    if run_btn:
        template = clean_seq(seq_input)
        if not template:
            st.error("Bitte eine DNA-Sequenz eingeben.")
            st.stop()

        args = {
            "PRIMER_TASK": "generic",
            "PRIMER_OPT_SIZE": primer_len,
            "PRIMER_MIN_SIZE": primer_len-2,
            "PRIMER_MAX_SIZE": primer_len+2,
            "PRIMER_OPT_TM": tm_opt,
            "PRIMER_MIN_TM": tm_opt-3,
            "PRIMER_MAX_TM": tm_opt+3,
            "PRIMER_MIN_GC": gc_min,
            "PRIMER_MAX_GC": gc_max,
            "PRIMER_SALT_MONOVALENT": salt_mM,
            "PRIMER_SALT_DIVALENT": mg_mM,
            "PRIMER_DNA_CONC": 50,
            "PRIMER_PRODUCT_SIZE_RANGE": [[int(a), int(b)] for a,b in (x.split("-") for x in product_size.split(",") if "-" in x)],
            "PRIMER_NUM_RETURN": 5
        }
        if pick_probe:
            args["PRIMER_PICK_INTERNAL_OLIGO"] = 1
            args["PRIMER_INTERNAL_OPT_TM"] = tm_opt + 8
            args["PRIMER_INTERNAL_MIN_TM"] = tm_opt + 5
            args["PRIMER_INTERNAL_MAX_TM"] = tm_opt + 10

        res = primer3.bindings.designPrimers({"SEQUENCE_ID":"template","SEQUENCE_TEMPLATE":template}, args)

        n_pairs = res.get("PRIMER_PAIR_NUM_RETURNED", 0)
        if n_pairs == 0:
            st.warning("Keine Primerpaare gefunden – Parameter anpassen.")
            st.stop()

        data = []
        for i in range(n_pairs):
            fwd = res[f"PRIMER_LEFT_{i}_SEQUENCE"]
            rev = res[f"PRIMER_RIGHT_{i}_SEQUENCE"]
            prod = res[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"]
            score, m = score_pair(fwd, rev)
            probe = res.get(f"PRIMER_INTERNAL_{i}_SEQUENCE") if pick_probe else None

            st.markdown(f"### 📌 Paar #{i+1} – Produkt {prod} bp · AI-Score {score}/100")
            st.code(f"Forward: {fwd}\nReverse: {rev}")
            st.pyplot(plot_dg(m))
            if probe:
                st.markdown(f"**Probe:** {reporter} — `{probe}` — {quencher}")

            data.append({
                "Pair": i+1,"Forward":fwd,"Reverse":rev,
                "Amplicon(bp)":prod,"Score":score,
                "Tm_F":m["tm_f"],"Tm_R":m["tm_r"],
                "ΔTm":m["tm_diff"],"GC_F":m["gc_f"],
                "GC_R":m["gc_r"],"ΔG_cross":m["dG_cross"],
                "Probe":probe or "-"
            })

        df = pd.DataFrame(data)
        st.subheader("📋 Zusammenfassung")
        st.dataframe(df, use_container_width=True)

        csv = df.to_csv(index=False).encode("utf-8")
        st.download_button("⬇️ CSV exportieren", csv, "primer_results.csv","text/csv")

        fasta_buf = StringIO()
        for _,r in df.iterrows():
            fasta_buf.write(f">Pair{r['Pair']}_F\n{r['Forward']}\n>Pair{r['Pair']}_R\n{r['Reverse']}\n")
            if r["Probe"] != "-":
                fasta_buf.write(f">Pair{r['Pair']}_Probe\n{r['Probe']}\n")
        st.download_button("⬇️ FASTA exportieren", fasta_buf.getvalue(), "primers.fasta","text/plain")

        st.success("Berechnung abgeschlossen ✅")
