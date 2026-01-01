
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D

def generate_plasmid_animation():
    print("Initializing HAWRA Plasmid 3D Animation...")
    
    # Create figure
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    fig.patch.set_facecolor('black')
    ax.set_facecolor('black')
    
    # Remove axes for cleaner look
    ax.set_axis_off()
    
    # --- Plasmid Torus ---
    theta = np.linspace(0, 2*np.pi, 100)
    phi = np.linspace(0, 2*np.pi, 50)
    theta, phi = np.meshgrid(theta, phi)
    
    R, r = 2, 0.6  # Major and minor radii
    
    # Initial coordinates
    x = (R + r*np.cos(phi)) * np.cos(theta)
    y = (R + r*np.cos(phi)) * np.sin(theta)
    z = r*np.sin(phi)
    
    # Plot plasmid backbone
    # Using plot_surface for a solid look or scatter for particle look
    # Let's use scatter for "digital/quantum" feel
    points_theta = np.linspace(0, 2*np.pi, 200)
    points_phi = np.linspace(0, 2*np.pi, 20)
    pt_theta, pt_phi = np.meshgrid(points_theta, points_phi)
    
    px = (R + r*np.cos(pt_phi)) * np.cos(pt_theta)
    py = (R + r*np.cos(pt_phi)) * np.sin(pt_theta)
    pz = r*np.sin(pt_phi)
    
    # Flatten arrays for scatter
    px = px.flatten()
    py = py.flatten()
    pz = pz.flatten()
    
    # Color mapping (gradient)
    colors = plt.cm.viridis(np.linspace(0, 1, len(px)))
    
    scatter = ax.scatter(px, py, pz, c=colors, s=2, alpha=0.6)
    
    # --- Functional Modules (Genes) ---
    # Define positions for modules on the torus
    modules = [
        {'name': 'P700 (Qubit)', 'angle': 0, 'color': 'red'},
        {'name': 'Lsi1 (Shield)', 'angle': np.pi/2, 'color': 'cyan'},
        {'name': 'PhyB (Input)', 'angle': np.pi, 'color': 'green'},
        {'name': 'Luc (Readout)', 'angle': 3*np.pi/2, 'color': 'orange'}
    ]
    
    module_scatters = []
    module_texts = []
    
    for mod in modules:
        angle = mod['angle']
        mx = (R + 1.2*r) * np.cos(angle)
        my = (R + 1.2*r) * np.sin(angle)
        mz = 0
        
        # Glow effect (large point)
        s = ax.scatter([mx], [my], [mz], c=mod['color'], s=200, alpha=0.8, edgecolors='white')
        module_scatters.append(s)
        
        # Label
        t = ax.text(mx, my, mz+0.5, mod['name'], color='white', fontsize=10, ha='center')
        module_texts.append(t)

    # Title
    ax.set_title("HAWRA Plasmid: Computational Operation", color='white', fontsize=16)
    
    # --- Animation Update ---
    def update(frame):
        # Rotate view
        ax.view_init(elev=30, azim=frame)
        
        # Pulsing effect for modules (simulate activity)
        # P700 pulses fast (quantum), others slower
        
        # Update colors/sizes randomly to simulate "processing"
        current_sizes = [200 + 50*np.sin(frame/10 + i) for i in range(len(modules))]
        
        # Update logic would go here if scatter returned a handle we could easily update sizes for
        # Matplotlib 3D scatter size update is tricky, so we'll just rotate for now
        # and maybe flash the title or add a dynamic text
        
        state_text = "Status: IDLE"
        if frame % 100 < 25:
            state_text = "Status: EXCITON TRANSFER (P700)"
        elif frame % 100 < 50:
            state_text = "Status: SILICA SHIELDING (Lsi1)"
        elif frame % 100 < 75:
            state_text = "Status: GENE REGULATION (Hill)"
        else:
            state_text = "Status: READOUT (Luciferase)"
            
        ax.set_title(f"HAWRA Metabiotic OS\n{state_text}", color='white', fontsize=14)
        
        return scatter,
    
    # Create animation
    print("Rendering frames (this may take a minute)...")
    anim = FuncAnimation(fig, update, frames=np.arange(0, 360, 2), interval=50)
    
    # Save
    output_path = 'results/hawra_plasmid_operation.gif'
    anim.save(output_path, writer='pillow', fps=20)
    print(f"Animation saved to {output_path}")

if __name__ == "__main__":
    generate_plasmid_animation()
