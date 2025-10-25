import streamlit as st
import io

def render():
    st.header("Reports")
    st.caption("Export quick reports as Markdown / Schneller Report-Export als Markdown")

    report_md = st.text_area("Report content (Markdown)", height=200,
                             value="# AI Primer Design Pro Report\n\n- Project: Demo\n- Date: TBD\n- Notes: ...")
    if st.button("Download .md"):
        bio = io.BytesIO(report_md.encode("utf-8"))
        st.download_button("Download report.md", data=bio.getvalue(), file_name="report.md", mime="text/markdown")