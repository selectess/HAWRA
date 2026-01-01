import matplotlib.pyplot as plt
from dna_features_viewer import GraphicFeature, CircularGraphicRecord
from Bio import SeqIO
import sys

# --- Charger le vrai plasmide depuis le fichier GenBank ---
genbank_file = "01_genomics/plasmids/HAWRA_FINAL_VALIDATED.gb"
record = SeqIO.read(genbank_file, "genbank")

# --- Créer un enregistrement graphique à partir des données réelles ---
graphic_record = CircularGraphicRecord(sequence_length=len(record.seq), features=record.features)

# --- Personnaliser les couleurs et les labels pour la clarté ---
def get_feature_color(feature):
    gene_colors = {
        "psaA": "#4b0082",
        "CRY2": "#0000ff",
        "Luc": "#00ff00",
        "Lsi1": "#ffA500",
        "HSP70": "#ff4500",
        "PEPC": "#dc143c",
        "P700": "#ff1493",
    }
    gene_name = feature.qualifiers.get('gene', [''])[0]
    if gene_name in gene_colors:
        return gene_colors[gene_name]
    # Couleur par défaut pour les autres features
    if feature.type == 'CDS':
        return 'cyan'
    return "grey"

def get_feature_label(feature):
    """Extract a single string label for the feature."""
    gene = feature.qualifiers.get('gene')
    if gene:
        return gene[0]
    locus_tag = feature.qualifiers.get('locus_tag')
    if locus_tag:
        return locus_tag[0]
    return feature.type

# Extraire TOUTES les features, pas seulement les CDS
all_features = [
    GraphicFeature(
        start=f.location.start, 
        end=f.location.end, 
        strand=f.location.strand, 
        color=get_feature_color(f),
        label=get_feature_label(f)
    )
    for f in record.features
]

# Vérification : s'assurer que des features ont été extraites
if not all_features:
    print("ERREUR : Aucune caractéristique (feature) n'a été trouvée dans le fichier GenBank.", file=sys.stderr)
    print("L'image générée sera probablement vide.", file=sys.stderr)
    sys.exit(1)

graphic_record.features = all_features

# --- Tracer le plasmide validé ---
ax, _ = graphic_record.plot(figure_width=10)
ax.set_title(f"Visualisation Validée du Plasmide HAWRA ({record.id}) - Toutes les features", weight="bold")
output_path = "05_data/results/hawra_plasmid_validated_visualization.png"
ax.figure.savefig(output_path, dpi=300)

print(f"Visualisation validée du plasmide HAWRA sauvegardée dans {output_path}")
print(f"{len(all_features)} caractéristiques ont été visualisées.")
