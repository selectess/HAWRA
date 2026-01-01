from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import os

# Définir les gènes et leurs fichiers FASTA
genes = {
    "35S": "V00111",
    "psaA": "NC_031161.1",
    "CRY2": "AF156319.1",
    "Luc": "U47132.1",
    "Lsi1": "KT159333.1",
    "HSP70": "AY842476.2",
    "PEPC": "X64138"
}

# Répertoire des séquences
input_dir = "01_genomics/raw_sequences/CDS"
output_dir = "02_synthetic_biology/constructs"
os.makedirs(output_dir, exist_ok=True)

# Charger les séquences
seqs = {}
for name, acc in genes.items():
    fasta_file = os.path.join(input_dir, f"{name}.fasta")
    if os.path.exists(fasta_file):
        seqs[name] = SeqIO.read(fasta_file, "fasta").seq
        print(f"Séquence {name} chargée : {len(seqs[name])} bp")
    else:
        print(f"Fichier pour {name} non trouvé : {fasta_file}")

# Mutation pour Luc rouge (S284T)
if "Luc" in seqs:
    luc_red = seqs["Luc"]
    # La mutation S284T est à la position 850-852 (codon TCT -> ACT)
    luc_red = luc_red[:850] + Seq("ACT") + luc_red[853:]
    print("Mutation S284T introduite dans le gène Luc.")
else:
    luc_red = Seq("")

# Ordre d'assemblage de la cassette
cassette_order = ["35S", "psaA", "CRY2", "Lsi1", "HSP70", "PEPC"]
final_cassette = Seq("")

for gene_name in cassette_order:
    if gene_name in seqs:
        final_cassette += seqs[gene_name]

# Ajouter le gène Luc muté et un terminateur simple
final_cassette += luc_red
final_cassette += Seq("TAG") # Terminateur simple

# Créer l'enregistrement GenBank
final_record = SeqRecord(
    final_cassette,
    id="HAWRA_PQPE_v1",
    name="HAWRA_Plasmid",
    description="Plasmide synthétique pour HAWRA avec qubits stables et géométrie optimisée"
)
final_record.annotations["molecule_type"] = "DNA"

# Ajouter des annotations de fonctionnalités (Features)
from Bio.SeqFeature import SeqFeature, FeatureLocation

features = [
    (seqs.get("35S", ""), "35S_promoter", "promoter"),
    (seqs.get("psaA", ""), "psaA_qubit_core", "CDS"),
    (seqs.get("CRY2", ""), "CRY2_quantum_gate", "CDS"),
    (luc_red, "Luc_red_reporter", "CDS"),
    (seqs.get("Lsi1", ""), "Lsi1_silica_cage", "CDS"),
    (seqs.get("HSP70", ""), "HSP70_thermal_protection", "CDS"),
    (seqs.get("PEPC", ""), "PEPC_energy_source", "CDS"),
]

current_pos = 0
for seq, name, ftype in features:
    if len(seq) > 0:
        location = FeatureLocation(start=current_pos, end=current_pos + len(seq))
        feature = SeqFeature(location, type=ftype, qualifiers={"label": name})
        final_record.features.append(feature)
        current_pos += len(seq)

# Sauvegarder le fichier GenBank
output_file = os.path.join(output_dir, "HAWRA_PLASMID.gb")
SeqIO.write(final_record, output_file, "genbank")

print(f"\nPlasmide HAWRA généré : {output_file}")
print(f"Longueur totale de la cassette : {len(final_cassette)} bp")
