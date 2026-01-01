#!/usr/bin/env python3
"""
HAWRA Demo Pipeline for Publication Figures
Generates publication-ready figures demonstrating bio-quantum computing capabilities
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import FancyBboxPatch, Circle, Rectangle
import json
import os
from datetime import datetime

class HAWRAFigures:
    """Generate publication figures for HAWRA bio-quantum computing framework"""
    
    def __init__(self):
        self.figures_dir = "publication_figures"
        os.makedirs(self.figures_dir, exist_ok=True)
        plt.style.use('seaborn-v0_8-whitegrid')
        
    def figure1_conceptual_overview(self):
        """Figure 1: Conceptual overview of HAWRA framework"""
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
        fig.suptitle('HAWRA Framework: PhytoQuantum Computing for Gene Regulatory Networks', 
                    fontsize=16, fontweight='bold')
        
        # Panel A: PQPE Architecture
        ax1.set_title('A. PhytoQuantum Processing Entity (PQPE)', fontsize=14, fontweight='bold')
        
        # Draw PQPE components
        components = [
            ('Quantum Core', (0.2, 0.8), '#FF6B6B'),
            ('Bio Interface', (0.8, 0.8), '#4ECDC4'),
            ('GRN Engine', (0.2, 0.2), '#45B7D1'),
            ('Light Controller', (0.8, 0.2), '#96CEB4')
        ]
        
        for name, (x, y), color in components:
            circle = Circle((x, y), 0.15, color=color, alpha=0.7)
            ax1.add_patch(circle)
            ax1.text(x, y, name, ha='center', va='center', fontsize=10, fontweight='bold')
        
        # Add connections
        for i, (name1, pos1, _) in enumerate(components):
            for name2, pos2, _ in components[i+1:]:
                ax1.plot([pos1[0], pos2[0]], [pos1[1], pos2[1]], 'k--', alpha=0.5)
        
        ax1.set_xlim(0, 1)
        ax1.set_ylim(0, 1)
        ax1.set_aspect('equal')
        ax1.axis('off')
        
        # Panel B: Arbol Language Syntax
        ax2.set_title('B. Arbol Language for Bio-Quantum Experiments', fontsize=14, fontweight='bold')
        
        arbol_code = """
# Gene Regulatory Network Definition
gene tf_a {
    promoter: strong
    rbs: medium
    cds: "TF_A"
    terminator: double
}

# Light Stimulus Coupling
stimulus light_pulse {
    wavelength: 660nm
    intensity: 100μmol/m²/s
    duration: 30min
}

# Quantum Operations
logical_qubit q1 is synthetic_biology_qubit
H on q1 with {
    coupling: tf_a.expression
    timescale: transcriptional
}

measure q1 -> result
        """
        
        ax2.text(0.05, 0.95, arbol_code, transform=ax2.transAxes, fontsize=10,
                verticalalignment='top', fontfamily='monospace',
                bbox=dict(boxstyle="round,pad=0.5", facecolor='#F8F9FA', alpha=0.8))
        ax2.set_xlim(0, 1)
        ax2.set_ylim(0, 1)
        ax2.axis('off')
        
        # Panel C: Hill Function Response Curves
        ax3.set_title('C. Hill Function Kinetics for Light-Responsive Genes', fontsize=14, fontweight='bold')
        
        # Hill function parameters
        K_values = [10, 50, 100]  # Different dissociation constants
        n_values = [2, 4, 6]      # Different Hill coefficients
        light_range = np.logspace(-1, 3, 100)
        
        colors = ['#FF6B6B', '#4ECDC4', '#45B7D1']
        
        for i, (K, n) in enumerate(zip(K_values, n_values)):
            response = (light_range**n) / (K**n + light_range**n)
            ax3.loglog(light_range, response, color=colors[i], linewidth=3,
                      label=f'K={K}, n={n}')
        
        ax3.set_xlabel('Light Intensity (μmol/m²/s)', fontsize=12)
        ax3.set_ylabel('Gene Expression Response', fontsize=12)
        ax3.legend(fontsize=10)
        ax3.grid(True, alpha=0.3)
        
        # Panel D: BSIM Compilation Pipeline
        ax4.set_title('D. BSIM Bytecode Compilation Pipeline', fontsize=14, fontweight='bold')
        
        pipeline_stages = [
            ('Arbol Source', '#FF6B6B'),
            ('Lexer', '#4ECDC4'),
            ('Parser', '#45B7D1'),
            ('Code Gen', '#96CEB4'),
            ('BSIM JSON', '#FECA57')
        ]
        
        y_pos = np.arange(len(pipeline_stages))
        
        for i, (stage, color) in enumerate(pipeline_stages):
            rect = Rectangle((0, i*0.8), 1, 0.6, color=color, alpha=0.7)
            ax4.add_patch(rect)
            ax4.text(0.5, i*0.8 + 0.3, stage, ha='center', va='center', 
                    fontsize=12, fontweight='bold')
            
            if i < len(pipeline_stages) - 1:
                ax4.arrow(0.5, i*0.8, 0, -0.2, head_width=0.05, head_length=0.05, 
                         fc='black', ec='black')
        
        ax4.set_xlim(-0.1, 1.1)
        ax4.set_ylim(-0.5, 4.5)
        ax4.axis('off')
        
        plt.tight_layout()
        plt.savefig(f'{self.figures_dir}/figure1_conceptual_overview.png', dpi=300, bbox_inches='tight')
        plt.close()
        
    def figure2_grn_simulation(self):
        """Figure 2: Gene Regulatory Network Simulation Results"""
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
        fig.suptitle('Gene Regulatory Network Simulation: Light-Responsive Circuit', 
                    fontsize=16, fontweight='bold')
        
        # Panel A: Network Topology
        ax1.set_title('A. GRN Topology with Light-Responsive Elements', fontsize=14, fontweight='bold')
        
        # Network nodes
        nodes = {
            'Light': (0.5, 0.9),
            'PhyB': (0.2, 0.7),
            'PIF3': (0.8, 0.7),
            'TF_A': (0.2, 0.4),
            'TF_B': (0.8, 0.4),
            'Output': (0.5, 0.1)
        }
        
        # Draw nodes
        for name, (x, y) in nodes.items():
            color = '#FF6B6B' if name == 'Light' else '#4ECDC4' if name in ['PhyB', 'PIF3'] else '#45B7D1'
            circle = Circle((x, y), 0.08, color=color, alpha=0.7)
            ax1.add_patch(circle)
            ax1.text(x, y, name, ha='center', va='center', fontsize=10, fontweight='bold')
        
        # Draw edges
        edges = [
            ('Light', 'PhyB'),
            ('Light', 'PIF3'),
            ('PhyB', 'TF_A'),
            ('PIF3', 'TF_B'),
            ('TF_A', 'Output'),
            ('TF_B', 'Output')
        ]
        
        for start, end in edges:
            x1, y1 = nodes[start]
            x2, y2 = nodes[end]
            ax1.arrow(x1, y1, x2-x1, y2-y1, head_width=0.03, head_length=0.03, 
                     fc='black', ec='black', alpha=0.6)
        
        ax1.set_xlim(0, 1)
        ax1.set_ylim(0, 1)
        ax1.set_aspect('equal')
        ax1.axis('off')
        
        # Panel B: Time-Dynamics Simulation
        ax2.set_title('B. Time-Dynamics of Light-Induced Gene Expression', fontsize=14, fontweight='bold')
        
        time = np.linspace(0, 24, 1000)  # 24 hours
        light_schedule = np.where((time % 12) < 6, 1, 0)  # 6h light, 6h dark cycles
        
        # Simulate gene expression dynamics
        def gene_expression(light, K=50, n=4, baseline=0.1):
            return baseline + (1 - baseline) * (light**n) / (K**n + light**n)
        
        tf_a_expr = gene_expression(light_schedule * 100)  # Light intensity 100 μmol/m²/s
        tf_b_expr = gene_expression(light_schedule * 80)  # Slightly different response
        
        ax2.plot(time, light_schedule * 100, 'k--', alpha=0.5, label='Light Intensity')
        ax2.plot(time, tf_a_expr, '#FF6B6B', linewidth=3, label='TF_A Expression')
        ax2.plot(time, tf_b_expr, '#4ECDC4', linewidth=3, label='TF_B Expression')
        
        ax2.set_xlabel('Time (hours)', fontsize=12)
        ax2.set_ylabel('Expression Level', fontsize=12)
        ax2.legend(fontsize=10)
        ax2.grid(True, alpha=0.3)
        
        # Panel C: Dose-Response Analysis
        ax3.set_title('C. Dose-Response Curves for Light-Responsive Genes', fontsize=14, fontweight='bold')
        
        light_doses = np.logspace(0, 3, 50)  # 1 to 1000 μmol/m²/s
        
        # Different Hill coefficients for different genes
        tf_a_response = gene_expression(light_doses, K=50, n=4)
        tf_b_response = gene_expression(light_doses, K=80, n=2)
        output_response = gene_expression(light_doses, K=60, n=3)
        
        ax3.semilogx(light_doses, tf_a_response, '#FF6B6B', linewidth=3, label='TF_A')
        ax3.semilogx(light_doses, tf_b_response, '#4ECDC4', linewidth=3, label='TF_B')
        ax3.semilogx(light_doses, output_response, '#45B7D1', linewidth=3, label='Output')
        
        ax3.set_xlabel('Light Intensity (μmol/m²/s)', fontsize=12)
        ax3.set_ylabel('Normalized Response', fontsize=12)
        ax3.legend(fontsize=10)
        ax3.grid(True, alpha=0.3)
        
        # Panel D: Sensitivity Analysis
        ax4.set_title('D. Sensitivity Analysis of GRN Parameters', fontsize=14, fontweight='bold')
        
        # Parameter ranges
        K_range = np.linspace(10, 200, 20)
        n_range = np.linspace(1, 8, 20)
        
        K_mesh, n_mesh = np.meshgrid(K_range, n_range)
        
        # Calculate sensitivity (derivative of response with respect to parameters)
        light_test = 100  # Test light intensity
        response = (light_test**n_mesh) / (K_mesh**n_mesh + light_test**n_mesh)
        
        sensitivity = np.abs(np.gradient(response, K_range, n_range)[0])
        
        contour = ax4.contourf(K_mesh, n_mesh, sensitivity, levels=20, cmap='viridis')
        ax4.set_xlabel('Dissociation Constant K', fontsize=12)
        ax4.set_ylabel('Hill Coefficient n', fontsize=12)
        
        cbar = plt.colorbar(contour, ax=ax4)
        cbar.set_label('Sensitivity', fontsize=12)
        
        plt.tight_layout()
        plt.savefig(f'{self.figures_dir}/figure2_grn_simulation.png', dpi=300, bbox_inches='tight')
        plt.close()
        
    def figure3_quantum_bio_hybrid(self):
        """Figure 3: Quantum-Bio Hybrid Simulation Results"""
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
        fig.suptitle('Quantum-Bio Hybrid Simulation: Light-Coupled Quantum Operations', 
                    fontsize=16, fontweight='bold')
        
        # Panel A: Quantum Circuit Diagram
        ax1.set_title('A. Quantum Circuit with Bio-Coupled Gates', fontsize=14, fontweight='bold')
        
        # Draw quantum circuit
        qubit_lines = [0.8, 0.6, 0.4, 0.2]
        
        # Qubit labels
        for i, label in enumerate(['|q0⟩', '|q1⟩', '|q2⟩', '|q3⟩']):
            ax1.text(-0.1, qubit_lines[i], label, fontsize=12, fontweight='bold')
        
        # Draw qubit lines
        for y in qubit_lines:
            ax1.plot([0, 1], [y, y], 'k', linewidth=2)
        
        # Gate positions
        gates = [
            ('H', 0.1, 0),  # Hadamard on q0
            ('X', 0.2, 1),  # Pauli-X on q1
            ('CNOT', 0.3, (0, 2)),  # CNOT from q0 to q2
            ('Bio-H', 0.5, 3),  # Bio-coupled Hadamard on q3
            ('M', 0.8, 0),  # Measurement on q0
            ('M', 0.8, 1),  # Measurement on q1
        ]
        
        for gate, x, qubit in gates:
            if isinstance(qubit, tuple):  # Two-qubit gate
                control, target = qubit
                # Control dot
                ax1.plot(x, qubit_lines[control], 'ko', markersize=8)
                # Target cross
                ax1.plot([x-0.03, x+0.03], [qubit_lines[target]-0.03, qubit_lines[target]+0.03], 'k-', linewidth=2)
                ax1.plot([x-0.03, x+0.03], [qubit_lines[target]+0.03, qubit_lines[target]-0.03], 'k-', linewidth=2)
                # Vertical line
                ax1.plot([x, x], [qubit_lines[control], qubit_lines[target]], 'k-', linewidth=2)
            else:
                # Single-qubit gate
                rect = Rectangle((x-0.04, qubit_lines[qubit]-0.04), 0.08, 0.08, 
                               facecolor='#FF6B6B' if 'Bio' in gate else '#4ECDC4', 
                               edgecolor='black', linewidth=2)
                ax1.add_patch(rect)
                ax1.text(x, qubit_lines[qubit], gate, ha='center', va='center', 
                        fontsize=10, fontweight='bold')
        
        ax1.set_xlim(-0.2, 1.2)
        ax1.set_ylim(0.1, 0.9)
        ax1.axis('off')
        
        # Panel B: Quantum State Evolution
        ax2.set_title('B. Quantum State Evolution with Bio-Coupling', fontsize=14, fontweight='bold')
        
        time = np.linspace(0, 10, 1000)  # 10 time units
        
        # Simulate quantum state evolution with bio-coupling
        def quantum_evolution(t, coupling_strength=0.1, frequency=1.0):
            # Rabi oscillations with bio-coupling modulation
            omega = frequency * (1 + coupling_strength * np.sin(0.5 * t))
            return np.cos(omega * t)**2, np.sin(omega * t)**2
        
        prob_0, prob_1 = quantum_evolution(time)
        
        ax2.plot(time, prob_0, '#FF6B6B', linewidth=3, label='|0⟩ Probability')
        ax2.plot(time, prob_1, '#4ECDC4', linewidth=3, label='|1⟩ Probability')
        
        # Add bio-coupling envelope
        bio_envelope = 0.1 * np.sin(0.5 * time)
        ax2.fill_between(time, 0, 1, alpha=0.2, color='#96CEB4', 
                        label='Bio-Coupling Modulation')
        
        ax2.set_xlabel('Time (arbitrary units)', fontsize=12)
        ax2.set_ylabel('State Probability', fontsize=12)
        ax2.legend(fontsize=10)
        ax2.grid(True, alpha=0.3)
        
        # Panel C: Classical Bit Correlation
        ax3.set_title('C. Classical Bit Correlation with Gene Expression', fontsize=14, fontweight='bold')
        
        # Simulate measurement outcomes
        n_measurements = 100
        measurement_times = np.linspace(0, 10, n_measurements)
        
        # Gene expression levels (affects measurement basis)
        gene_expr = 0.5 + 0.3 * np.sin(measurement_times)
        
        # Quantum measurement outcomes influenced by gene expression
        random_outcomes = np.random.random(n_measurements)
        quantum_outcomes = (random_outcomes < 0.5 + 0.2 * gene_expr).astype(int)
        
        # Plot correlation
        ax3.scatter(gene_expr, quantum_outcomes, c=measurement_times, 
                   cmap='viridis', alpha=0.7, s=50)
        
        # Add trend line
        z = np.polyfit(gene_expr, quantum_outcomes, 1)
        p = np.poly1d(z)
        ax3.plot(gene_expr, p(gene_expr), "r--", alpha=0.8, linewidth=2)
        
        ax3.set_xlabel('Gene Expression Level', fontsize=12)
        ax3.set_ylabel('Quantum Measurement Outcome', fontsize=12)
        ax3.set_ylim(-0.1, 1.1)
        
        cbar = plt.colorbar(ax3.collections[0], ax=ax3)
        cbar.set_label('Time', fontsize=12)
        
        # Panel D: Bio-Quantum Fidelity
        ax4.set_title('D. Bio-Quantum Operation Fidelity', fontsize=14, fontweight='bold')
        
        # Simulate fidelity over time with different coupling strengths
        coupling_strengths = [0.05, 0.1, 0.2, 0.3]
        colors = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#96CEB4']
        
        time_fidelity = np.linspace(0, 20, 100)
        
        for i, coupling in enumerate(coupling_strengths):
            # Fidelity decays with time but improves with coupling
            baseline_fidelity = 0.95
            decay_rate = 0.01 * (1 + coupling)
            coupling_benefit = coupling * (1 - np.exp(-time_fidelity / 5))
            
            fidelity = baseline_fidelity * np.exp(-decay_rate * time_fidelity) + coupling_benefit
            fidelity = np.clip(fidelity, 0, 1)  # Ensure fidelity stays in [0,1]
            
            ax4.plot(time_fidelity, fidelity, color=colors[i], linewidth=3,
                    label=f'Coupling = {coupling}')
        
        ax4.set_xlabel('Time (arbitrary units)', fontsize=12)
        ax4.set_ylabel('Operation Fidelity', fontsize=12)
        ax4.legend(fontsize=10)
        ax4.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(f'{self.figures_dir}/figure3_quantum_bio_hybrid.png', dpi=300, bbox_inches='tight')
        plt.close()
        
    def figure4_3d_molecular(self):
        """Figure 4: 3D Molecular Visualization of Bio-Quantum Interactions"""
        fig = plt.figure(figsize=(20, 16))
        
        # Create a 2x2 grid with different 3D views
        gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.3)
        
        fig.suptitle('3D Molecular Visualization: Bio-Quantum Interface', 
                    fontsize=20, fontweight='bold')
        
        # Panel A: Plasmid DNA Structure
        ax1 = fig.add_subplot(gs[0, 0], projection='3d')
        ax1.set_title('A. Plasmid DNA with Light-Responsive Promoters', fontsize=16, fontweight='bold')
        
        # Generate toroidal plasmid structure
        theta = np.linspace(0, 2*np.pi, 100)
        phi = np.linspace(0, 2*np.pi, 50)
        theta, phi = np.meshgrid(theta, phi)
        
        # Torus parameters
        R, r = 2, 0.8  # Major and minor radii
        
        x = (R + r*np.cos(phi)) * np.cos(theta)
        y = (R + r*np.cos(phi)) * np.sin(theta)
        z = r*np.sin(phi)
        
        # Color by base pair position
        colors = plt.cm.viridis(np.linspace(0, 1, len(theta)))
        
        for i in range(0, len(theta), 5):
            ax1.plot(x[i], y[i], z[i], 'o', color=colors[i], markersize=4, alpha=0.7)
        
        # Mark light-responsive regions
        light_regions = [(0, 20), (40, 60), (80, 100)]
        for start, end in light_regions:
            ax1.plot(x[start:end], y[start:end], z[start:end], 'ro', 
                    markersize=8, alpha=0.8, label='Light-Responsive' if start == 0 else "")
        
        ax1.set_xlabel('X (nm)', fontsize=12)
        ax1.set_ylabel('Y (nm)', fontsize=12)
        ax1.set_zlabel('Z (nm)', fontsize=12)
        ax1.legend(fontsize=10)
        
        # Panel B: Protein-DNA Interactions
        ax2 = fig.add_subplot(gs[0, 1], projection='3d')
        ax2.set_title('B. Transcription Factor Binding to DNA', fontsize=16, fontweight='bold')
        
        # DNA double helix
        t = np.linspace(0, 4*np.pi, 200)
        helix_radius = 1
        dna_x = t
        dna_y1 = helix_radius * np.cos(t)
        dna_z1 = helix_radius * np.sin(t)
        dna_y2 = helix_radius * np.cos(t + np.pi)
        dna_z2 = helix_radius * np.sin(t + np.pi)
        
        ax2.plot(dna_x, dna_y1, dna_z1, 'b-', linewidth=4, alpha=0.8, label='DNA Strand 1')
        ax2.plot(dna_x, dna_y2, dna_z2, 'g-', linewidth=4, alpha=0.8, label='DNA Strand 2')
        
        # Transcription factors
        tf_positions = [np.pi, 3*np.pi, 5*np.pi, 7*np.pi]
        for i, pos in enumerate(tf_positions):
            tf_x = pos
            tf_y = 2 * np.cos(pos)
            tf_z = 2 * np.sin(pos)
            
            # Draw TF as sphere
            u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
            tf_sphere_x = tf_x + 0.5 * np.cos(u) * np.sin(v)
            tf_sphere_y = tf_y + 0.5 * np.sin(u) * np.sin(v)
            tf_sphere_z = tf_z + 0.5 * np.cos(v)
            
            ax2.plot_surface(tf_sphere_x, tf_sphere_y, tf_sphere_z, 
                           color='red', alpha=0.6)
        
        ax2.set_xlabel('Helix Axis', fontsize=12)
        ax2.set_ylabel('Y (nm)', fontsize=12)
        ax2.set_zlabel('Z (nm)', fontsize=12)
        ax2.legend(fontsize=10)
        
        # Panel C: Quantum Dot-Plasmid Coupling
        ax3 = fig.add_subplot(gs[1, 0], projection='3d')
        ax3.set_title('C. Quantum Dot Coupled to Plasmid', fontsize=16, fontweight='bold')
        
        # Quantum dot
        u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
        qd_radius = 0.5
        qd_x = qd_radius * np.cos(u) * np.sin(v)
        qd_y = qd_radius * np.sin(u) * np.sin(v)
        qd_z = qd_radius * np.cos(v)
        
        ax3.plot_surface(qd_x, qd_y, qd_z, color='purple', alpha=0.8, label='Quantum Dot')
        
        # Plasmid (smaller torus)
        theta_small = np.linspace(0, 2*np.pi, 50)
        phi_small = np.linspace(0, 2*np.pi, 25)
        theta_small, phi_small = np.meshgrid(theta_small, phi_small)
        
        R_small, r_small = 1, 0.3
        
        plasmid_x = 2 + (R_small + r_small*np.cos(phi_small)) * np.cos(theta_small)
        plasmid_y = (R_small + r_small*np.cos(phi_small)) * np.sin(theta_small)
        plasmid_z = r_small*np.sin(phi_small)
        
        ax3.plot_surface(plasmid_x, plasmid_y, plasmid_z, 
                       color='orange', alpha=0.6, label='Plasmid')
        
        # Coupling field
        coupling_field = np.sqrt(qd_x**2 + (qd_y-1.5)**2 + qd_z**2)
        ax3.contour(qd_x, qd_y, qd_z, levels=[0.8], colors='cyan', linewidths=2)
        
        ax3.set_xlabel('X (nm)', fontsize=12)
        ax3.set_ylabel('Y (nm)', fontsize=12)
        ax3.set_zlabel('Z (nm)', fontsize=12)
        ax3.legend(fontsize=10)
        
        # Panel D: Base Pairing Interactions
        ax4 = fig.add_subplot(gs[1, 1], projection='3d')
        ax4.set_title('D. Watson-Crick Base Pairing', fontsize=16, fontweight='bold')
        
        # Base pairs
        base_pairs = [
            ('A', 'T', 0, 'red', 'blue'),
            ('T', 'A', 1, 'blue', 'red'),
            ('G', 'C', 2, 'green', 'orange'),
            ('C', 'G', 3, 'orange', 'green'),
        ]
        
        for i, (base1, base2, pos, color1, color2) in enumerate(base_pairs):
            z = pos * 0.5
            
            # Base 1
            ax4.scatter(i, 0, z, c=color1, s=200, alpha=0.8)
            ax4.text(i, 0, z, base1, ha='center', va='center', fontsize=12, fontweight='bold')
            
            # Base 2
            ax4.scatter(i, 1, z, c=color2, s=200, alpha=0.8)
            ax4.text(i, 1, z, base2, ha='center', va='center', fontsize=12, fontweight='bold')
            
            # Hydrogen bonds
            ax4.plot([i, i], [0, 1], [z, z], 'k--', linewidth=2, alpha=0.6)
            
            # Base stacking interactions
            if i > 0:
                ax4.plot([i-1, i], [0, 0], [z-0.5, z], 'k:', linewidth=1, alpha=0.4)
                ax4.plot([i-1, i], [1, 1], [z-0.5, z], 'k:', linewidth=1, alpha=0.4)
        
        ax4.set_xlabel('Base Pair Position', fontsize=12)
        ax4.set_ylabel('Strand', fontsize=12)
        ax4.set_zlabel('Helix Axis', fontsize=12)
        
        plt.tight_layout()
        plt.savefig(f'{self.figures_dir}/figure4_3d_molecular.png', dpi=300, bbox_inches='tight')
        plt.close()
        
    def generate_all_figures(self):
        """Generate all publication figures"""
        print("Generating Figure 1: Conceptual Overview...")
        self.figure1_conceptual_overview()
        
        print("Generating Figure 2: GRN Simulation...")
        self.figure2_grn_simulation()
        
        print("Generating Figure 3: Quantum-Bio Hybrid...")
        self.figure3_quantum_bio_hybrid()
        
        print("Generating Figure 4: 3D Molecular Visualization...")
        self.figure4_3d_molecular()
        
        print(f"All figures generated in '{self.figures_dir}' directory")
        
        # Create figure summary
        summary = {
            "figures": [
                {
                    "name": "figure1_conceptual_overview.png",
                    "description": "Conceptual overview of HAWRA framework showing PQPE architecture, Arbol language syntax, Hill function kinetics, and BSIM compilation pipeline",
                    "panels": ["PQPE Architecture", "Arbol Language", "Hill Functions", "BSIM Pipeline"]
                },
                {
                    "name": "figure2_grn_simulation.png", 
                    "description": "Gene Regulatory Network simulation results showing light-responsive circuit topology, time-dynamics, dose-response curves, and sensitivity analysis",
                    "panels": ["GRN Topology", "Time-Dynamics", "Dose-Response", "Sensitivity Analysis"]
                },
                {
                    "name": "figure3_quantum_bio_hybrid.png",
                    "description": "Quantum-Bio hybrid simulation results showing quantum circuits with bio-coupled gates, state evolution, classical bit correlation, and operation fidelity",
                    "panels": ["Quantum Circuit", "State Evolution", "Bit Correlation", "Operation Fidelity"]
                },
                {
                    "name": "figure4_3d_molecular.png",
                    "description": "3D molecular visualization showing plasmid DNA structure, protein-DNA interactions, quantum dot coupling, and Watson-Crick base pairing",
                    "panels": ["Plasmid Structure", "Protein-DNA Binding", "Quantum Dot Coupling", "Base Pairing"]
                }
            ],
            "generated_at": datetime.now().isoformat(),
            "total_figures": 4
        }
        
        with open(f'{self.figures_dir}/figures_summary.json', 'w') as f:
            json.dump(summary, f, indent=2)
        
        return summary

def main():
    """Main function to generate all publication figures"""
    print("HAWRA Demo Pipeline for Publication Figures")
    print("=" * 50)
    
    figures = HAWRAFigures()
    summary = figures.generate_all_figures()
    
    print("\nPublication figures generated successfully!")
    print(f"Total figures: {summary['total_figures']}")
    print(f"Output directory: {figures.figures_dir}")
    
    # Print figure details
    for fig_info in summary['figures']:
        print(f"\n- {fig_info['name']}: {fig_info['description'][:80]}...")

if __name__ == "__main__":
    main()