
from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np
import os

def visualize_plasmid(gb_file, output_image):
    """
    Génère une carte circulaire du plasmide à partir d'un fichier GenBank.
    """
    record = SeqIO.read(gb_file, "genbank")
    plasmid_len = len(record)

    fig, ax = plt.subplots(figsize=(12, 12), subplot_kw=dict(projection="polar"))
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    ax.spines['polar'].set_visible(False)

    # Gènes et leurs couleurs
    gene_colors = {
        "psaA": "red",
        "CRY2": "blue",
        "SIT1": "green",
        "Lsi1": "green",
        "LUC": "orange",
        "HSP70": "purple",
        "PEPC": "brown"
    }

    # Dessiner les features utiles (promoter, 5'UTR, CDS, terminator)
    for feature in record.features:
        ftype = feature.type
        start = int(feature.location.start)
        end = int(feature.location.end)
        theta_start = (start / plasmid_len) * 2 * np.pi
        theta_end = (end / plasmid_len) * 2 * np.pi

        label = None
        color = "gray"
        bottom = 0.8
        height = 0.1

        if ftype == "promoter":
            label = feature.qualifiers.get("label", ["promoter"])[0]
            color = "gold"
            bottom = 0.65
            height = 0.08
        elif ftype == "misc_feature":
            qlabel = feature.qualifiers.get("label", ["misc"])[0]
            label = qlabel
            color = "teal" if "5UTR" in qlabel else "cyan"
            bottom = 0.72
            height = 0.06
        elif ftype == "CDS":
            gene_name = feature.qualifiers.get("gene", ["unknown"])[0]
            label = gene_name
            gn_upper = gene_name.upper()
            color = gene_colors.get(gene_name, gene_colors.get(gn_upper, "gray"))
            bottom = 0.8
            height = 0.12
        elif ftype == "terminator":
            label = feature.qualifiers.get("label", ["terminator"])[0]
            color = "darkred"
            bottom = 0.6
            height = 0.06
        else:
            continue

        ax.bar(x=(theta_start + theta_end) / 2,
               height=height,
               width=theta_end - theta_start,
               bottom=bottom,
               color=color,
               alpha=0.7)

        angle = ((theta_start + theta_end) / 2) * 180 / np.pi
        if 90 < angle < 270:
            angle -= 180
        ax.text(x=(theta_start + theta_end) / 2,
                y=bottom + height + 0.02,
                s=label,
                rotation=angle,
                ha='center',
                va='center',
                fontsize=9,
                fontweight='bold')

    # Ajouter le titre
    ax.text(0, 0, f"Carte du Plasmide HAWRA\n{plasmid_len} bp", 
            ha='center', va='center', fontsize=16, fontweight='bold')

# Créer le répertoire de sortie s'il n'existe pas
    output_dir = os.path.dirname(output_image)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    plt.savefig(output_image, dpi=300, bbox_inches='tight')
    print(f"Carte du plasmide sauvegardée dans {output_image}")

if __name__ == "__main__":
    visualize_plasmid(
        "05_data/results/hawra_plasmid.gb",
        "05_data/results/genetic_visualization/hawra_plasmid_map.png"
    )
