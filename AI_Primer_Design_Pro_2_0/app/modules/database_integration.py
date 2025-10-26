# ==========================================
# 🧫 Database & Reference Integration Module
# ==========================================
import streamlit as st
import requests
import sqlite3
import os

# Lokaler SQLite Cache (wird automatisch erstellt)
DB_PATH = os.path.join(os.path.dirname(__file__), "local_cache.db")

def init_cache():
    conn = sqlite3.connect(DB_PATH)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS api_cache (
            source TEXT,
            query TEXT,
            result TEXT,
            PRIMARY KEY (source, query)
        )
    """)
    conn.commit()
    conn.close()

# ---------- Funktionen ----------
def query_ncbi(term, db="gene"):
    """NCBI Entrez-Abfrage"""
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {"db": db, "term": term, "retmode": "json"}
    r = requests.get(url, params=params)
    return r.json() if r.status_code == 200 else {"error": "No result"}

def query_uniprot(uniprot_id):
    """UniProt-Daten abrufen"""
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.txt"
    r = requests.get(url)
    return r.text if r.status_code == 200 else "Kein Eintrag gefunden."

def query_neb_enzyme(enzyme):
    """NEB Restriktionsenzyme (API-ähnlich)"""
    fake_db = {
        "EcoRI": {"Recognition": "GAATTC", "Cleavage": "G^AATTC", "Source": "E. coli"},
        "BamHI": {"Recognition": "GGATCC", "Cleavage": "G^GATCC", "Source": "B. amyloliquefaciens"},
        "HindIII": {"Recognition": "AAGCTT", "Cleavage": "A^AGCTT", "Source": "H. influenzae"},
    }
    return fake_db.get(enzyme, {"Error": "Enzym nicht gefunden."})

# ---------- Streamlit Interface ----------
def run_database_integration():
    st.title("🧫 Database & Reference Integration")
    st.caption("NCBI · UniProt · NEB Enzyme · SQLite Cache · Cloud Sync")

    init_cache()

    tabs = st.tabs([
        "🔬 NCBI Entrez",
        "🧬 UniProt",
        "🧫 NEB Restriction Enzymes",
        "💾 Local SQLite Cache",
        "☁️ Cloud Sync (API Key)",
    ])

    # --- NCBI ---
    with tabs[0]:
        st.subheader("🔬 NCBI Entrez-Abfragen")
        db = st.selectbox("Datenbank wählen", ["gene", "protein", "nucleotide"])
        term = st.text_input("Suchbegriff (z. B. BRCA1, TP53)")
        if st.button("Suchen", key="ncbi"):
            res = query_ncbi(term, db)
            st.json(res)

    # --- UniProt ---
    with tabs[1]:
        st.subheader("🧬 UniProt Integration")
        uid = st.text_input("UniProt ID (z. B. P69905)")
        if st.button("Abrufen", key="uniprot"):
            data = query_uniprot(uid)
            st.text_area("Ergebnis:", data, height=300)

    # --- NEB Enzyme ---
    with tabs[2]:
        st.subheader("🧫 NEB Restriktionsenzyme")
        enzyme = st.text_input("Enzymname (z. B. EcoRI)")
        if st.button("Suchen", key="neb"):
            info = query_neb_enzyme(enzyme)
            st.json(info)

    # --- Local Cache ---
    with tabs[3]:
        st.subheader("💾 Lokaler SQLite Cache")
        if st.button("Cache anzeigen"):
            conn = sqlite3.connect(DB_PATH)
            df = None
            try:
                df = st.dataframe(conn.execute("SELECT * FROM api_cache").fetchall())
            except Exception:
                st.info("Cache ist leer.")
            conn.close()

    # --- Cloud Sync ---
    with tabs[4]:
        st.subheader("☁️ Cloud-Sync mit API-Key Login")
        api_key = st.text_input("API-Key", type="password")
        if api_key:
            st.success("API-Key gespeichert (Demo-Modus).")
        else:
            st.info("Bitte API-Key eingeben, um Cloud-Funktion zu aktivieren.")
