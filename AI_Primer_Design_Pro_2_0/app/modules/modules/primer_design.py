# modules/primer_design.py

import streamlit as st
import primer3
from Bio.Seq import Seq
from Bio.Data.IUPACData import ambiguous_dna_values
import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt
from io import StringIO

# ----------------------------------------
# 🔧 FUNKTION: Modulstart
# ----------------------------------------
def run():
    st.title("🧬 Primer Design & PCR Tools")
    st.caption("Primer3 integriert • Dimer/Hairpin/Self/Cross • Degenerate • ΔG-Analyse • qPCR-Probes • Annotation & Amplicon Preview")

    # ---------- Hilfsfunktionen ----------
    def clean_seq(s: str) -> str:
        """Bereinigt FASTA/Plain-Sequenzen"""
        s = s.strip()
        if s.startswith(">"):
            s = "\n".join([line for line in s.splitlines() if not line.startswith(">")])
        s = re.sub(r"[^ACGTURYKMSWBDHVNacgturykmswbdhvn]", "", s)
        return s.upper()

    def annotate_sequence(seq, f_start, f_len, r_start, r_len, amplicon_start, amplicon_len):
        """ASCII-Minimap mit Primer/amplicon-Markern."""
        n = len(seq)
        line = ["·"] * n
        for i in range(f_start, min(f_start + f_len, n)):
            line[i] = "F"
        for i in range(r_start, min(r_start + r_len, n)):
            line[i] = "R"
        for i in range(amplicon_start, min(amplicon_start + amplicon_len, n)):
            if line[i] == "·":
                line[i] = "-"
        return "".join(line)

    def dimer_hairpin_scores(fwd, rev):
        sd_f = primer3.calcHomodimer(fwd).dg
        hp_f = primer3.calcHairpin(fwd).dg
        sd_r = primer3.calcHomodimer(rev).dg
        hp_r = primer3.calcHairpin(rev).dg
        cd = primer3.calcHeterodimer(fwd, rev).dg
        return sd_f, hp_f, sd_r, hp_r, cd

    def simple_efficiency_score(fwd, rev):
        """Heuristischer Score (0–100)"""
        def feat(p):
            gc = (p.count("G")+p.count("C"))/len(p)*100
            tm = primer3.calcTm(p, mv_conc=50, dv_conc=1.5/2, dntp_conc=0.8, dna_conc=50)
            hp = primer3.calcHairpin(p).dg
            sd = primer3.calcHomodimer(p).dg
            return gc, tm, hp, sd

        gcf, tmf, hpf, sdf = feat(fwd)
        gcr, tmr, hpr, sdr = feat(rev)
        cd = primer3.calcHeterodimer(fwd, rev).dg

        gc_ok = 40 <= (gcf+gcr)/2 <= 60
        tm_diff = abs(tmf - tmr)
        tm_ok = tm_diff <= 2.5
        hp_pen = max(0, (-min(hpf, hpr) - 3)) * 5
        sd_pen = max(0, (-min(sdf, sdr) - 5)) * 3
        cd_pen = max(0, (-cd - 6)) * 4

        base = 85
        if not gc_ok: base -= 10
        if not tm_ok: base -= 10
        score = base - hp_pen - sd_pen - cd_pen
        return int(np.clip(score, 0, 100)), {
            "tm_f": round(tmf,1), "tm_r": round(tmr,1),
            "tm_diff": round(tm_diff,1), "cd_dG": round(cd,1),
            "hp_f": round(hpf,1), "hp_r": round(hpr,1),
            "sd_f": round(sdf,1), "sd_r": round(sdr,1),
            "gc_f": round(gcf,1), "gc_r": round(gcr,1)
        }

    def plot_dg_bars(metrics):
        fig, ax = plt.subplots()
        names = ["hp_f","sd_f","hp_r","sd_r","cross"]
        vals = [metrics["hp_f"], metrics["sd_f"], metrics["hp_r"], metrics["sd_r"], metrics["cd_dG"]]
        ax.bar(names, vals)
        ax.set_ylabel("ΔG (kcal/mol)")
        ax.set_title("Thermodynamische Bewertung")
        return fig

    # ---------- Sidebar ----------
    with st.sidebar:
        st.header("⚙️ Parameter")
        seq_input = st.text_area("Template-Sequenz (FASTA oder Text)", height=150)
        product_size = st.text_input("Produktgröße (z. B. 80–200)", value="80-200")
        primer_len = st.slider("Primerlänge (optimal)", 18, 30, 20)
        tm_opt = st.slider("Tm optimal (°C)", 55, 72, 60)
        gc_min = st.number_input("GC% min", 20, 80, 40)
        gc_max = st.number_input("GC% max", 20, 80, 60)
        salt_mM = st.number_input("Salz [Na⁺] (mM)", 1.0, 200.0, 50.0)
        mg_mM = st.number_input("Mg²⁺ (mM)", 0.0, 5.0, 1.5)
        pick_internal = st.checkbox("qPCR-Probe (TaqMan / MB) hinzufügen", value=True)
        reporter = st.selectbox("Reporter", ["FAM","HEX","ROX","Cy5","VIC"])
        quencher = st.selectbox("Quencher", ["BHQ1","BHQ2","TAMRA"])
        run_btn = st.button("🚀 Primer entwerfen")

    # ---------- Hauptbereich ----------
    if run_btn:
        template = clean_seq(seq_input)
        if not template:
            st.error("Bitte eine DNA-Template-Sequenz eingeben.")
            st.stop()

        with st.spinner("Primer3 wird ausgeführt…"):
            args = {
                'PRIMER_TASK': 'generic',
                'PRIMER_OPT_SIZE': primer_len,
                'PRIMER_MIN_SIZE': primer_len-2,
                'PRIMER_MAX_SIZE': primer_len+2,
                'PRIMER_OPT_TM': tm_opt,
                'PRIMER_MIN_TM': tm_opt-3,
                'PRIMER_MAX_TM': tm_opt+3,
                'PRIMER_MIN_GC': gc_min,
                'PRIMER_MAX_GC': gc_max,
                'PRIMER_SALT_MONOVALENT': salt_mM,
                'PRIMER_SALT_DIVALENT': mg_mM,
                'PRIMER_DNA_CONC': 50,
                'PRIMER_PRODUCT_SIZE_RANGE': [[int(x.split('-')[0]), int(x.split('-')[1])] for x in product_size.split(',') if '-' in x],
                'PRIMER_NUM_RETURN': 5
            }

            if pick_internal:
                args['PRIMER_PICK_INTERNAL_OLIGO'] = 1
                args['PRIMER_INTERNAL_OPT_TM'] = tm_opt + 8
                args['PRIMER_INTERNAL_MIN_TM'] = tm_opt + 5
                args['PRIMER_INTERNAL_MAX_TM'] = tm_opt + 10

            res = primer3.bindings.designPrimers(
                {'SEQUENCE_ID': 'template', 'SEQUENCE_TEMPLATE': template},
                args
            )

        n_pairs = res.get('PRIMER_PAIR_NUM_RETURNED', 0)
        if n_pairs == 0:
            st.warning("Keine Primerpaare gefunden – bitte Parameter anpassen.")
            st.stop()

        data = []
        for i in range(n_pairs):
            fwd = res[f'PRIMER_LEFT_{i}_SEQUENCE']
            rev = res[f'PRIMER_RIGHT_{i}_SEQUENCE']
            f_pos, f_len = res[f'PRIMER_LEFT_{i}']
            r_pos, r_len = res[f'PRIMER_RIGHT_{i}']
            prod = res[f'PRIMER_PAIR_{i}_PRODUCT_SIZE']

            sd_f, hp_f, sd_r, hp_r, cd = dimer_hairpin_scores(fwd, rev)
            score, metr = simple_efficiency_score(fwd, rev)
            amp_line = annotate_sequence(template, f_pos, f_len, r_pos, r_len, f_pos, prod)

            st.markdown(f"### 📌 Primer-Paar #{i+1}  (Amplicon: {prod} bp, AI-Score: {score}/100)")
            st.code(f"Forward: {fwd}\nReverse: {rev}", language="plaintext")
            st.text_area("Amplicon-Preview", amp_line, height=100)
            st.pyplot(plot_dg_bars(metr))

            probe_seq = res.get(f'PRIMER_INTERNAL_{i}_SEQUENCE') if pick_internal else None
            if probe_seq:
                st.markdown(f"**Probe:** {reporter} — `{probe_seq}` — {quencher}")

            data.append({
                "Pair": i+1, "Forward": fwd, "Reverse": rev,
                "Amplicon(bp)": prod, "AI-Score": score,
                "Tm_F": metr["tm_f"], "Tm_R": metr["tm_r"],
                "ΔG_cross": metr["cd_dG"], "GC_F": metr["gc_f"], "GC_R": metr["gc_r"]
            })

        df = pd.DataFrame(data)
        st.markdown("### 📋 Zusammenfassung aller Primer")
        st.dataframe(df, use_container_width=True)

        # Export
        csv = df.to_csv(index=False).encode("utf-8")
        st.download_button("⬇️ CSV exportieren", csv, "primer_design_results.csv", "text/csv")

        fasta_buf = StringIO()
        for _, r in df.iterrows():
            fasta_buf.write(f">Primer_{r['Pair']}_F\n{r['Forward']}\n>Primer_{r['Pair']}_R\n{r['Reverse']}\n")
        st.download_button("⬇️ FASTA exportieren", fasta_buf.getvalue(), "primers.fasta", "text/plain")

        st.success("Analyse abgeschlossen ✅")

