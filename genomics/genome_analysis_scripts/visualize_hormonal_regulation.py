
import matplotlib.pyplot as plt
import os

def visualize_hormonal_regulation(output_image):
    """
    Génère une visualisation conceptuelle de la régulation hormonale dans le projet HAWRA.
    """
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 7)
    ax.axis('off')

    # Titre
    fig.suptitle('Régulation Hormonale Conceptuelle pour HAWRA', fontsize=16, fontweight='bold')

    # Nœuds
    nodes = {
        'Lumière': (1, 5),
        'Horloge Circadienne': (3, 5),
        'Production d\'Auxine': (5, 5),
        'Croissance de la Plante': (7, 5),
        'Production de Qubits': (9, 5),
        'Stress Hydrique': (1, 2),
        'Production d\'ABA': (3, 2),
        'Fermeture des Stomates': (5, 2),
        'Conservation de l\'Eau': (7, 2),
    }

    # Dessiner les nœuds
    for node, (x, y) in nodes.items():
        ax.text(x, y, node, ha='center', va='center', 
                bbox=dict(boxstyle='round,pad=0.5', fc='lightblue', alpha=0.7))

    # Flèches
    arrows = [
        ('Lumière', 'Horloge Circadienne'),
        ('Horloge Circadienne', 'Production d\'Auxine'),
        ('Production d\'Auxine', 'Croissance de la Plante'),
        ('Croissance de la Plante', 'Production de Qubits'),
        ('Stress Hydrique', 'Production d\'ABA'),
        ('Production d\'ABA', 'Fermeture des Stomates'),
        ('Fermeture des Stomates', 'Conservation de l\'Eau'),
    ]

    for start_node, end_node in arrows:
        start_pos = nodes[start_node]
        end_pos = nodes[end_node]
        ax.annotate('', xy=end_pos, xytext=start_pos, 
                    arrowprops=dict(arrowstyle='->', lw=2, color='gray'))

    # Créer le répertoire de sortie s'il n'existe pas
    output_dir = os.path.dirname(output_image)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    plt.savefig(output_image, dpi=300, bbox_inches='tight')
    print(f"Graphique de régulation hormonale sauvegardé dans {output_image}")

if __name__ == "__main__":
    visualize_hormonal_regulation(
        "/Users/mehdiwhb/Desktop/HAWRA/05_data/results/hormonal_regulation/hormonal_regulation_concept.png"
    )
