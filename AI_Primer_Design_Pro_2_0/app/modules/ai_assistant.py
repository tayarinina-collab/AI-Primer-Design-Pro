import openai
import os

def interpret_sequence(seq, lang="DE"):
    prompt_de = f"""
    Du bist ein KI-Laborassistent. Analysiere die folgende Sequenz und erkläre sie in einfachem Labor-Deutsch:
    {seq}
    Gib eine kurze, klare Zusammenfassung ihrer biologischen Bedeutung.
    """

    prompt_en = f"""
    You are an AI lab assistant. Analyze the following sequence and explain it in simple laboratory English:
    {seq}
    """

    try:
        openai.api_key = os.getenv("OPENAI_API_KEY")
        response = openai.ChatCompletion.create(
            model="gpt-4o-mini",
            messages=[{"role": "system", "content": prompt_de if lang == "DE" else prompt_en}]
        )
        return response["choices"][0]["message"]["content"]
    except Exception:
        # Offline fallback
        return "⚠️ Offline-Modus aktiv – KI-Analyse momentan nicht verfügbar."
