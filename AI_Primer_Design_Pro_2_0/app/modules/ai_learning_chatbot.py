# -*- coding: utf-8 -*-
"""
AI Learning & Chatbot System
- AI Lab Assistant Chatbot
- Adaptive Learning
- Explainable AI Reports
- Knowledge Base Mode
"""

import streamlit as st
import random

def run_ai_learning_chatbot():
    st.title("🤖 AI Learning & Chatbot System")
    st.caption("Intelligenter Laborassistent mit Lernfähigkeit und Erklärfunktionen")

    # -----------------------------
    # Session Memory (Lernmodus)
    # -----------------------------
    if "user_prefs" not in st.session_state:
        st.session_state["user_prefs"] = {
            "preferred_tm": 60,
            "gc_range": (40, 60),
            "language": "Deutsch"
        }

    st.markdown("### 🧠 Adaptive Learning")
    st.write("Aktuell gespeicherte Benutzerpräferenzen:")
    st.json(st.session_state["user_prefs"])

    if st.button("🔁 Präferenzen zurücksetzen"):
        st.session_state["user_prefs"] = {
            "preferred_tm": 60,
            "gc_range": (40, 60),
            "language": "Deutsch"
        }
        st.success("✅ Zurückgesetzt!")

    st.markdown("---")

    # -----------------------------
    # Chatbot Interface
    # -----------------------------
    st.subheader("💬 AI Lab Assistant Chatbot")
    st.caption("Stelle eine Frage zu Primerdesign, PCR, oder Laborprotokollen:")

    user_input = st.text_area("🧪 Deine Frage:", placeholder="Warum funktioniert mein PCR-Ansatz nicht?")
    mode = st.radio(
        "Chatbot-Modus wählen:",
        ["Lab Assistant", "Explainable AI", "Knowledge Base"],
        horizontal=True
    )

    if st.button("🚀 Antwort generieren"):
        if not user_input.strip():
            st.warning("Bitte eine Frage eingeben.")
        else:
            response = generate_ai_response(user_input, mode)
            st.markdown("### 🧬 Antwort")
            st.write(response)

            # Simulierte Lernfunktion: Nutzerpräferenzen anpassen
            if "Tm" in user_input or "Temperatur" in user_input:
                st.session_state["user_prefs"]["preferred_tm"] = random.choice(range(57, 64))
            if "GC" in user_input or "Gehalt" in user_input:
                st.session_state["user_prefs"]["gc_range"] = (random.randint(35, 45), random.randint(55, 65))

    st.markdown("---")
    st.caption("🧠 AI Primer Design Pro · Adaptive Chat System · Version 1.0")


# -----------------------------
# Antwortgenerator (Demo)
# -----------------------------
def generate_ai_response(question, mode):
    """Offline-Dummy-Chatbot mit themenspezifischen Antworten"""
    q = question.lower()

    # Lab Assistant Mode
    if mode == "Lab Assistant":
        if "primer" in q or "design" in q:
            return (
                "🔬 **Design-Tipp:** Überprüfe die GC-Verteilung. "
                "Ein stabiler Primer liegt meist zwischen 40–60 % GC. "
                "Achte außerdem auf keine Hairpins oder Dimerbildung."
            )
        elif "pcr" in q:
            return (
                "🧪 **PCR-Fehleranalyse:** Prüfe die Annealing-Temperatur. "
                "Falls keine Banden sichtbar sind, versuche eine niedrigere Temperatur (−2 °C). "
                "Kontrolliere auch Mg²⁺-Konzentration und Enzymaktivität."
            )
        else:
            return (
                "🤖 Ich bin dein Laborassistent! "
                "Frag mich alles zu PCR, Primerdesign oder Klonierung – ich erkläre jeden Schritt."
            )

    # Explainable AI Mode
    elif mode == "Explainable AI":
        return (
            "📊 **Erklärung:** Der Primer wurde als *suboptimal* eingestuft, "
            "weil sein ΔG-Wert unter −6 kcal/mol liegt (Hinweis auf sekundäre Strukturen). "
            "Außerdem ist der GC-Gehalt außerhalb des bevorzugten Bereichs (40–60 %)."
        )

    # Knowledge Base Mode
    elif mode == "Knowledge Base":
        if "taq" in q:
            return "🧬 **Taq-Polymerase**: Ein thermostabiles Enzym aus *Thermus aquaticus*, essentiell für PCR."
        elif "annealing" in q:
            return "🌡️ **Annealing:** Phase der PCR, in der Primer an die DNA binden (typisch 55–65 °C)."
        elif "dntp" in q:
            return "💡 **dNTPs:** Bausteine der DNA-Synthese – Adenin, Cytosin, Guanin und Thymin-Nukleotide."
        else:
            return "📚 Ich kann wissenschaftliche Begriffe, Protokolle oder Enzyme erklären. Frag mich etwas Spezifisches!"

    return "⚙️ Keine Antwort gefunden – bitte versuche es mit einem anderen Thema."
