from Bio import Entrez, SeqIO
import os
import ssl

# Contournement de la vérification SSL
if not os.environ.get("PYTHONHTTPSVERIFY", "") and getattr(ssl, "_create_unverified_context", None):
    ssl._create_default_https_context = ssl._create_unverified_context

Entrez.email = "gemini@google.com"

genes = {
    "35S": "V00111",        # Promoteur constitutif
    "psaA": "NC_031161.1",  # Qubit cœur (cohérence quantique)
    "CRY2": "AF156319.1",   # Portes quantiques
    "Luc": "U47132.1",      # Lecture photonique
    "Lsi1": "KT159333.1",   # Cage de silice
    "HSP70": "AY842476.2",  # Protection thermique (nouvel identifiant)
    "PEPC": "X64138"        # Énergie 24h via CAM
}
# Créer le répertoire de sortie s'il n'existe pas
output_dir = "01_genomics/raw_sequences/CDS"
os.makedirs(output_dir, exist_ok=True)

for name, acc in genes.items():
    file_path = os.path.join(output_dir, f"{name}.fasta")
    if not os.path.exists(file_path):
        print(f"Téléchargement de {name} ({acc})...")
        try:
            handle = Entrez.efetch(db="nucleotide", id=acc, rettype="fasta", retmode="text")
            record = SeqIO.read(handle, "fasta")
            SeqIO.write(record, file_path, "fasta")
            print(f"{name} téléchargé : longueur {len(record.seq)} bp")
        except Exception as e:
            print(f"Erreur lors du téléchargement de {name}: {e}")
    else:
        record = SeqIO.read(file_path, "fasta")
        print(f"{name} existe déjà : longueur {len(record.seq)} bp")
