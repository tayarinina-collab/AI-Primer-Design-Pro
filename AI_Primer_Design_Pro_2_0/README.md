# AI Primer Design Pro (v2.0)

EN · An intelligent bioinformatics platform for sequence analysis, primer design, and molecular visualization.
DE · Eine intelligente Bioinformatik-App fuer Sequenzanalyse, Primer-Design und molekulare Visualisierung.

---

## English

AI Primer Design Pro combines advanced primer design, in-silico PCR, plasmid visualization, and protein analysis with an integrated AI assistant.
Developed as a modern alternative to Geneious Prime, optimized for efficiency, automation, and bilingual accessibility (EN/DE).

### Key Features (v2.0 Beta)
- Primer Design (Primer3 integration)
- In-Silico PCR (virtual amplification)
- Plasmid Map (basic visualization)
- Restriction Tools (Bio.Restriction)
- Protein Tools + AI Chatbot (requires OPENAI_API_KEY)
- Bilingual UI (English/German)
- Cloud-ready (Streamlit Cloud)

> Roadmap: NCBI / UniProt API, 3D Structure Viewer, MSA/Phylogeny, Workflow Automation, Reports (PDF), Versioning.

### Quickstart (Local)
    pip install -r requirements.txt
    cd app
    streamlit run main.py

### Deploy to Streamlit Cloud
1) Push this folder to a public GitHub repo.
2) Go to https://share.streamlit.io -> Deploy an app -> select your repo.
3) Set secrets (optional):
   - OPENAI_API_KEY for AI chatbot
4) The app builds and gives you a public URL.

---

## Deutsch

AI Primer Design Pro vereint modernes Primer-Design, In-Silico-PCR, Plasmid-Visualisierung und Protein-Analyse mit einem integrierten KI-Assistenten.
Entwickelt als moderne Alternative zu Geneious Prime - mit Fokus auf Effizienz, Automatisierung und zweisprachiger Oberflaeche (DE/EN).

### Hauptfunktionen (v2.0 Beta)
- Primer-Design (Primer3-Integration)
- In-Silico-PCR (virtuelle Amplifikation)
- Plasmid-Karte (Basis-Visualisierung)
- Restriktions-Tools (Bio.Restriction)
- Protein-Tools + KI-Chatbot (benoetigt OPENAI_API_KEY)
- Zweisprachige UI (Deutsch/Englisch)
- Cloud-faehig (Streamlit Cloud)

### Schnellstart (Lokal)
    pip install -r requirements.txt
    cd app
    streamlit run main.py

### Deployment auf Streamlit Cloud
1) Diesen Ordner in ein oeffentliches GitHub-Repo pushen.
2) https://share.streamlit.io -> Deploy an app -> Repo auswaehlen.
3) Secrets setzen (optional):
   - OPENAI_API_KEY fuer KI-Chatbot
4) Streamlit erstellt automatisch die oeffentliche URL.