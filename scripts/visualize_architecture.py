import matplotlib.pyplot as plt
import matplotlib.patches as patches

def create_hawra_architecture_visualization():
    """Crée une visualisation architecturale du projet HAWRA."""
    fig, ax = plt.subplots(figsize=(16, 10))
    ax.set_title('Architecture Conceptuelle du Projet HAWRA', fontsize=20, fontweight='bold')
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 70)
    ax.axis('off')

    # Styles
    box_style = dict(boxstyle='round,pad=0.5', fc='lightblue', ec='steelblue', lw=2)
    arrow_style = dict(arrowstyle='-|>', color='gray', lw=2)
    component_style = dict(boxstyle='round,pad=0.5', fc='ivory', ec='darkkhaki', lw=2)

    # --- Composants Principaux ---
    # BioOS
    ax.text(50, 60, 'BioOS / Core System', ha='center', va='center', fontsize=16, bbox=box_style)

    # ARBOL Language
    ax.text(15, 45, 'Langage ARBOL\n(Interface Utilisateur)', ha='center', va='center', fontsize=12, bbox=component_style)
    ax.add_patch(patches.FancyArrowPatch((25, 45), (40, 58), connectionstyle="arc3,rad=.3", **arrow_style))

    # Plasmide HAWRA (Wetware)
    ax.text(50, 30, 'Wetware: Plasmide HAWRA', ha='center', va='center', fontsize=14, bbox=dict(boxstyle='round,pad=0.5', fc='lightgreen', ec='darkgreen', lw=2))
    ax.add_patch(patches.FancyArrowPatch((50, 58), (50, 33), connectionstyle="arc3,rad=0", **arrow_style))

    # Bio-Qubit P700
    ax.text(50, 10, 'Unité de Calcul: Bio-Qubit P700', ha='center', va='center', fontsize=12, bbox=component_style)
    ax.add_patch(patches.FancyArrowPatch((50, 27), (50, 13), connectionstyle="arc3,rad=0", **arrow_style))

    # Simulations (Validation)
    ax.text(85, 45, 'Validation Numérique', ha='center', va='center', fontsize=14, bbox=dict(boxstyle='round,pad=0.5', fc='mistyrose', ec='salmon', lw=2))
    ax.add_patch(patches.FancyArrowPatch((60, 58), (75, 47), connectionstyle="arc3,rad=-.3", **arrow_style))

    # --- Sous-composants de validation ---
    ax.text(85, 35, 'Sim 1: Décohérence P700', ha='center', va='center', fontsize=10, bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='salmon'))
    ax.text(85, 25, 'Sim 2: Protocole PQPE', ha='center', va='center', fontsize=10, bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='salmon'))
    ax.add_patch(patches.FancyArrowPatch((85, 42), (85, 38), connectionstyle="arc3,rad=0", **arrow_style))
    ax.add_patch(patches.FancyArrowPatch((85, 32), (85, 28), connectionstyle="arc3,rad=0", **arrow_style))

    # Sauvegarde de l'image
    output_path = '/Users/mehdiwhb/Desktop/HAWRA/05_data/results/hawra_architecture.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Visualisation de l'architecture sauvegardée dans {output_path}")

if __name__ == "__main__":
    create_hawra_architecture_visualization()
