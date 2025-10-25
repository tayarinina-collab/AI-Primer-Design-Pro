import streamlit as st
try:
    import primer3
except Exception:
    primer3 = None

def render():
    st.header("Primer-Design")
    st.caption("Primer3 integration with sensible defaults / Primer3-Integration mit sinnvollen Defaults")

    seq = st.text_area("DNA sequence (5'->3') / DNA-Sequenz", height=150, value="ATGCGTACGTAGCTAGCTAGCTAGCTAGCATCGATCGATGCTAGCTAGCTAGC")
    target_start = st.number_input("Target start (bp) / Zielstart", min_value=0, value=10, step=1)
    target_len = st.number_input("Target length (bp) / Ziellaenge", min_value=1, value=120, step=1)

    col1, col2, col3 = st.columns(3)
    with col1:
        tm = st.slider("Tm (C)", 50, 70, (58, 62))
    with col2:
        gc = st.slider("GC %", 30, 80, (40, 60))
    with col3:
        plen = st.slider("Primer length", 16, 30, (20, 24))

    if st.button("Design primers / Primer entwerfen", use_container_width=True):
        if primer3 is None:
            st.error("primer3-py not available. Please install dependencies. (requirements.txt)")
            return

        seq_template = str(seq).upper()
        target = [int(target_start), int(target_len)]
        params = {
            "PRIMER_NUM_RETURN": 5,
            "PRIMER_OPT_SIZE": int(sum(plen)//2),
            "PRIMER_MIN_SIZE": int(plen[0]),
            "PRIMER_MAX_SIZE": int(plen[1]),
            "PRIMER_OPT_TM": sum(tm)/2,
            "PRIMER_MIN_TM": tm[0],
            "PRIMER_MAX_TM": tm[1],
            "PRIMER_MIN_GC": gc[0],
            "PRIMER_MAX_GC": gc[1],
            "PRIMER_PRODUCT_SIZE_RANGE": [[int(target_len*0.8), int(target_len*1.5)]],
        }

        res = primer3.bindings.designPrimers(
            {
                "SEQUENCE_ID": "target",
                "SEQUENCE_TEMPLATE": seq_template,
                "SEQUENCE_TARGET": target,
            },
            params
        )

        n = res.get("PRIMER_PAIR_NUM_RETURNED", 0)
        if n == 0:
            st.warning("No primers found with current constraints.")
        else:
            for i in range(n):
                left = res.get(f"PRIMER_LEFT_{i}_SEQUENCE", "")
                right = res.get(f"PRIMER_RIGHT_{i}_SEQUENCE", "")
                prod = res.get(f"PRIMER_PAIR_{i}_PRODUCT_SIZE", "?")
                st.write(f"Pair {i+1} - product {prod} bp")
                st.code(f"F: {left}\nR: {right}")
                import streamlit as st

def run_primer_design():
    st.title("🧬 Primer Design Modul")
    st.markdown("Dies ist das Primer Design Modul.")
