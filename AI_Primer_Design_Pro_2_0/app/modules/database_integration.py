# ==========================================
# 🧫 Database & Reference Integration Module
# ==========================================
import streamlit as st
import requests
import sqlite3
import os
import json
import pandas as pd

# ---------- Lokaler SQLite Cache ----------
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

# ---------- API-Funktionen ----------
def query_ncbi(term, db="gene"):
    """NCBI Entrez-Abfrage"""
    if not term.strip():
        return {"error": "Kein Suchbegriff eingegeben."}
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {"db": db, "term": term, "retmode": "json"}
    r = requests.get(url, params=params)
    if r.status_code == 200:
        return r.json()
    return {"error": "Keine Verbindung zur NCBI API."}

def query_uniprot(uniprot_id):
    """UniProt-Daten abrufen"""
    if not uniprot_id.strip():
        return "⚠️ Bitte eine UniProt-ID eingeben."
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.txt"
    r = requests.get(url)
    return r.text if r.status_code == 200 else "❌ Kein Eintrag gefunden."

def query_neb_enzyme(enzyme):
    """NEB Restriktionsenzyme (Demo-Datenbank)"""
    fake_db = {
        "EcoRI":  {"Recognition": "GAATTC", "Cleavage": "G^AATTC", "Source": "E. coli"},
        "BamHI":  {"Recognition": "GGATCC", "Cleavage": "G^GATCC", "Source": "B. amyloliquefaciens"},
        "HindIII":{"Recognition": "AAGCTT", "Cleavage": "A^AGCTT", "Source": "H. influenzae"},
    }
    return fake_db.get(enzyme.strip(), {"Error": "❌ Enzym nicht gefunden."})

# ---------- Hauptfunktion ----------
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

    # --- 🔬 NCBI ENTrez ---
    with tabs[0]:
        st.subheader("🔬 NCBI Entrez-Abfragen")

        db = st.selectbox("Datenbank wählen", ["gene", "protein", "nucleotide"])
        term = st.text_input("Suchbegriff (z. B. BRCA1, TP53)")

        if st.button("Suchen", key="ncbi"):
            res = query_ncbi(term, db)
            if "error" in res:
                st.error(res["error"])
            else:
                st.success("✅ Abfrage erfolgreich!")
                # schön formatiert
                st.json(res)
                # Kurzfassung anzeigen
                count = res.get("esearchresult", {}).get("count", "0")
                ids = res.get("esearchresult", {}).get("idlist", [])
                st.markdown(f"**📊 Treffer gesamt:** {count}")
                if ids:
                    st.markdown(f"**🧬 Erste IDs:** {', '.join(ids[:10])}")

    # --- 🧬 UniProt ---
    with tabs[1]:
        st.subheader("🧬 UniProt Integration")
        uid = st.text_input("UniProt ID (z. B. P69905)")

        if st.button("Abrufen", key="uniprot"):
            data = query_uniprot(uid)
            if data.startswith("⚠️") or data.startswith("❌"):
                st.warning(data)
            else:
                st.success("✅ UniProt-Daten erfolgreich geladen!")
                st.text_area("📄 Ergebnis (Rohdaten):", data, height=300)

    # --- 🧫 NEB Restriktionsenzyme ---
    with tabs[2]:
        st.subheader("🧫 NEB Restriktionsenzyme")
        enzyme = st.text_input("Enzymname (z. B. EcoRI)")

        if st.button("Suchen", key="neb"):
            info = query_neb_enzyme(enzyme)

            if "Error" in info:
                st.error(info["Error"])
            else:
                st.success(f"✅ Daten für {enzyme} gefunden:")
                # 1️⃣ JSON-Anzeige
                st.json(info)
                # 2️⃣ oder schöne Tabelle
                df = pd.DataFrame(info.items(), columns=["Eigenschaft", "Wert"])
                st.markdown("### 🧬 Enzym-Informationen")
                st.table(df)
                # 3️⃣ Kurztext mit Emojis
                st.markdown(f"""
                **Recognition site:** 🧩 `{info['Recognition']}`  
                **Cleavage pattern:** ✂️ `{info['Cleavage']}`  
                **Source organism:** 🧫 *{info['Source']}*
                """)

    # --- 💾 Lokaler Cache ---
    with tabs[3]:
        st.subheader("💾 Lokaler SQLite Cache")
        if st.button("Cache anzeigen"):
            conn = sqlite3.connect(DB_PATH)
            try:
                data = conn.execute("SELECT * FROM api_cache").fetchall()
                if data:
                    df = pd.DataFrame(data, columns=["Source", "Query", "Result"])
                    st.dataframe(df, use_container_width=True)
                else:
                    st.info("ℹ️ Cache ist leer.")
            except Exception as e:
                st.error(f"Fehler beim Lesen des Caches: {e}")
            conn.close()

    # --- ☁️ Cloud Sync ---
    with tabs[4]:
        st.subheader("☁️ Cloud-Sync mit API-Key Login")
        api_key = st.text_input("API-Key", type="password")
        if api_key:
            st.success("🔐 API-Key gespeichert (Demo-Modus).")
        else:
            st.info("Bitte API-Key eingeben, um Cloud-Funktion zu aktivieren.")
