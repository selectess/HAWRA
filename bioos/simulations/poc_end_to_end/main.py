import numpy as np
import qutip as qt
from qutip import Options
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import plotly.graph_objects as go

# --- 1. Modélisation du Bio-Qubit (P700) ---
def initialize_bio_qubit():
    """
    Initialise le bio-qubit dans son état fondamental |0>.
    Le P700 est modélisé comme un système à deux niveaux.
    """
    return qt.basis(2, 0)

# --- 2. Simulation de la Porte Hadamard (Impulsion Laser) ---
def apply_hadamard_gate(qubit_state):
    """
    Applique une porte de Hadamard au qubit pour le mettre en superposition.
    Ceci simule une impulsion laser π/2.
    """
    # Construction manuelle de la porte de Hadamard pour compatibilité
    h_matrix = (1 / np.sqrt(2)) * np.array([[1, 1], [1, -1]])
    h_gate = qt.Qobj(h_matrix, dims=[[2], [2]])
    return h_gate * qubit_state

# --- 3. Simulation de la Décohérence ---
def simulate_decoherence(initial_state, times, t1, t2):
    """
    Simule l'évolution du qubit soumis à la décohérence (T1 et T2).
    Utilise le solver `mesolve` de QuTiP.
    
    Args:
        initial_state: L'état quantique initial.
        times: Un tableau NumPy des points temporels pour la simulation.
        t1: Le temps de relaxation (amplitude damping).
        t2: Le temps de déphasage (phase damping).
    
    Returns:
        Le résultat de la simulation (objet `Result` de QuTiP).
    """
    # Opérateurs de collapse pour T1 et T2
    c_ops = []
    # T1 relaxation
    if t1 > 0:
        c_ops.append(np.sqrt(1/t1) * qt.sigmam())
    # T2 dephasing
    if t2 > 0:
        gamma_phi = 1/t2 - 1/(2*t1)
        if gamma_phi > 0:
            c_ops.append(np.sqrt(gamma_phi) * qt.sigmaz())
            
    # Hamiltonien (ici, nul car on observe la relaxation libre)
    H = 0 * qt.sigmaz()

    # Options pour forcer la sauvegarde des états
    options = Options(store_states=True)
    
    # Résolution de l'équation maîtresse
    result = qt.mesolve(H, initial_state, times, c_ops, e_ops=[qt.sigmax(), qt.sigmay(), qt.sigmaz()], options=options)
    return result

# --- 4. Lecture de l'État Final (Readout) ---
def readout(final_state):
    """
    Effectue une mesure dans la base Z.
    Traduit la probabilité d'être dans l'état |1> en une intensité lumineuse.
    """
    # Probabilité de mesurer l'état |1>
    prob_1 = qt.expect(qt.projection(2, 1, 1), final_state)
    
    # Modèle simple de conversion en intensité lumineuse (normalisée)
    # Plus la probabilité de |1> est élevée, plus la luminescence est forte
    light_intensity = prob_1
    return light_intensity

# --- 5. Visualisation ---
def visualize_simulation(result, times):
    """
    Génère deux visualisations :
    1. Une animation de la sphère de Bloch montrant la trajectoire du qubit.
    2. Un graphique de l'intensité lumineuse (readout) en fonction du temps.
    """
    # a) Animation de la sphère de Bloch
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    sphere = qt.Bloch(axes=ax)

    def animate(i):
        ax.clear()
        sphere.axes = ax
        sphere.clear()
        sphere.add_points([result.expect[0][:i+1], result.expect[1][:i+1], result.expect[2][:i+1]])
        sphere.make_sphere()
        return (ax,)

    ani = FuncAnimation(fig, animate, frames=len(times), blit=False)
    ani.save('bloch_sphere_decoherence.gif', writer='imagemagick', fps=20)
    plt.close(fig)
    print("Animation de la sphère de Bloch sauvegardée dans 'bloch_sphere_decoherence.gif'")

    # b) Graphique de l'intensité lumineuse en fonction du temps
    readout_over_time = [readout(state) for state in result.states]
    
    fig_intensity = go.Figure()
    fig_intensity.add_trace(go.Scatter(x=times, y=readout_over_time, mode='lines', name='Intensité Lumineuse Simulée'))
    fig_intensity.update_layout(
        title="Simulation du Readout par Bioluminescence",
        xaxis_title="Temps (ns)",
        yaxis_title="Intensité Lumineuse (normalisée)",
        template="plotly_dark"
    )
    fig_intensity.write_html("readout_intensity.html")
    print("Graphique de l'intensité lumineuse sauvegardé dans 'readout_intensity.html'")


# --- Main Orchestrator ---
def main():
    """
    Orchestre la simulation de bout en bout.
    """
    print("--- Début de la Simulation PoC HAWRA ---")
    
    # Paramètres de la simulation
    T1_decoherence_time = 1e-9   # Temps de relaxation en ns (supposé long)
    T2_decoherence_time = 650e-15 # Temps de déphasage en s (650 fs)
    simulation_time = 2e-12     # Durée totale de la simulation en s (2 ps)
    time_steps = 200            # Nombre de points temporels
    
    times = np.linspace(0, simulation_time, time_steps)
    
    # Étape 1: Initialisation
    print("1. Initialisation du bio-qubit P700...")
    qubit_initial = initialize_bio_qubit()
    
    # Étape 2: Application de la porte quantique
    print("2. Application de la porte Hadamard (superposition)...")
    qubit_superposition = apply_hadamard_gate(qubit_initial)
    
    # Étape 3: Simulation de la décohérence
    print(f"3. Simulation de la décohérence (T1={T1_decoherence_time*1e9:.1f}ns, T2={T2_decoherence_time*1e15:.1f}fs)...")
    result = simulate_decoherence(qubit_superposition, times, T1_decoherence_time, T2_decoherence_time)
    
    # Étape 4: Lecture (Readout)
    final_state = result.states[-1]
    intensity = readout(final_state)
    print(f"4. Lecture de l'état final : Intensité lumineuse = {intensity:.4f}")
    
    # Étape 5: Visualisation
    print("5. Génération des visualisations...")
    visualize_simulation(result, times)
    
    print("--- Simulation Terminée ---")
    print("Les résultats ont été sauvegardés dans 'bloch_sphere_decoherence.gif' et 'readout_intensity.html'.")

if __name__ == "__main__":
    main()