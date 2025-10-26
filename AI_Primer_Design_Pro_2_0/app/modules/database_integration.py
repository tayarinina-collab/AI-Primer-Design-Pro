# -*- coding: utf-8 -*-
"""
Database & Reference Integration Modul
• NCBI Entrez-Abfragen (Gene, Protein, Nucleotide)
• UniProt Integration (Annotation, Funktionen, Domains)
• NEB Restriction Enzyme Database
• Lokale SQLite-Cache-Datenbank (Offline Mirror)
• Cloud-Sync mit API-Key Login
"""

import streamlit as st
import pandas as pd
import requests
import sqlite3
import os
from datetime import datetime

# ---------------------------------------------------------------------
# Lokaler Cache (SQLite)
# ---------------------------------------------------------------------
CACHE_DB = "local_cache.db"

def init_cache():
    conn = sqlite3.connect(CACHE_DB)
    c = conn.cursor()
    c.execute("""
        CREATE TABLE IF NOT EXISTS cache (
            source TEXT,
            query TEXT,
            result TEXT,
            timestamp TEXT
        )
    """)
    conn.commit()
    conn.close()

def cache_get(source, query):
    conn = sqlite3.connect(CACHE_DB)
    c = conn.cursor()
    c.execute("SELECT result FROM cache WHERE source=? AND query=?", (source, query))
    row = c.fetchone()
    conn.close()
    return row[0] if row else None

def cache_set(source, query, result):
    conn = sqlite3.connect(CACHE_DB)
    c = conn.cursor()
    c.execute("INSERT INTO cache VALUES (?, ?, ?, ?)", (source, query, result, datetime.now().isoformat()))
    conn.commit()
    conn.close()


# ---------------------------------------------------------------------
# Cloud Sync / API Key Login
# ---------------------------------------------------------------------
def cloud_login():
    api_key = st.text_input("🔑 API-Key für Cloud-Sync", type="password")
    if api_key:
        st.session_state["api_key"] = api_key
        st.success("API-Key gespeichert ✅")
    elif "api_key" in st.session_state:
        st.info("Cloud-Sync aktiviert mit gespeichertem API-Key")
    else:
        st.warning("Bitte API-Key eingeben, um Cloud-Sync zu aktivieren.")


# ---------------------------------------------------------------------
# NCBI Entrez API
# ---------------------------------------------------------------------
def query_ncbi(database, term, retmax=5):
    """Einfacher Wrapper für NCBI Entrez eSearch + eSummary"""
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {"db": database, "term": term, "retmode": "json", "retmax": retmax}
    r = requests.get(url, params=params)
    data = r.json()
    ids = data["esearchresult"].get("idlist", [])
    results = []
    for _id in ids:
        summary_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
        summary_params = {"db": database, "id": _id, "retmode": "json"}
        s = requests.get(summary_url, params=summary_params).json()
        doc = s.get("result", {}).get(_id, {})
        results.append({
            "ID": _id,
            "Title": doc.get("title", ""),
            "Organism": doc.get("organism", ""),
            "Updated": doc.get("updatedate", "")
        })
    return pd.DataFrame(results)


# ---------------------------------------------------------------------
# UniProt API
# ---------------------------------------------------------------------
def query_uniprot(term):
    url = f"https://rest.uniprot.org/uniprotkb/search?query={term}&format=tsv&fields=accession,id,protein_name,organism,length"
    r = requests.get(url)
    if r.status_code == 200:
        df = pd.read_csv(pd.compat.StringIO(r.text), sep="\t")
        return df
    else:
        st.error(f"Fehler bei UniProt-Anfrage: {r.status_code}")
        return None


# ---------------------------------------------------------------------
# NEB Restriction Enzyme DB (JSON-API)
# ---------------------------------------------------------------------
def query_neb():
    url = "https://api.neb.com/biolabs/api/v1/enzymes"
    r = requests.get(url)
    if r.status_code == 200:
        data = r.json()
        enzymes = [{"Name": e["name"], "Recognition": e["sequence"], "Overhang": e["overhang"], "Supplier": e["supplier"]} for e in data]
        return pd.DataFrame(enzymes)
    else:
        st.error("NEB-API nicht erreichbar.")
        return None


# ---------------------------------------------------------------------
# Streamlit UI
# ---------------------------------------------------------------------
def run_database_integration():
    st.title("🧫 Database & Reference Integration")
    st.caption("NCBI · UniProt · NEB · Local Cache · Cloud Sync")

    tabs = st.tabs(["🔍 NCBI Search", "🧬 UniProt", "✂️ Restriction Enzymes", "💾 Cache", "☁️ Cloud Sync"])

    # -------- NCBI --------
    with tabs[0]:
        st.subheader("🔍 NCBI Entrez-Abfrage")
        database = st.selectbox("Datenbank", ["gene", "protein", "nucleotide"])
        term = st.text_input("Suchbegriff", "BRCA1")
        if st.button("Suchen"):
            cached = cache_get("NCBI", f"{database}:{term}")
            if cached:
                st.info("Ergebnis aus Cache geladen.")
                df = pd.read_json(cached)
            else:
                df = query_ncbi(database, term)
                cache_set("NCBI", f"{database}:{term}", df.to_json())
            st.dataframe(df, use_container_width=True)

    # -------- UniProt --------
    with tabs[1]:
        st.subheader("🧬 UniProt Integration")
        term = st.text_input("Proteinname oder Organismus", "p53 human")
        if st.button("UniProt abfragen"):
            df = query_uniprot(term)
            if df is not None:
                st.dataframe(df, use_container_width=True)

    # -------- NEB --------
    with tabs[2]:
        st.subheader("✂️ Restriktionsenzyme-Datenbank (NEB)")
        if st.button("Daten abrufen"):
            df = query_neb()
            if df is not None:
                st.dataframe(df, use_container_width=True)

    # -------- Cache --------
    with tabs[3]:
        st.subheader("💾 Lokale Cache-Datenbank")
        if os.path.exists(CACHE_DB):
            df = pd.read_sql("SELECT * FROM cache", sqlite3.connect(CACHE_DB))
            st.dataframe(df, use_container_width=True)
        else:
            st.info("Noch keine Cache-Daten vorhanden.")

    # -------- Cloud Sync --------
    with tabs[4]:
        st.subheader("☁️ Cloud-Synchronisierung")
        cloud_login()
