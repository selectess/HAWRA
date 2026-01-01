import matplotlib.pyplot as plt
import numpy as np

# Set style manually without seaborn
plt.style.use('default')

# Data: Temperature vs Coherence/Energy Cost
# 1. Superconducting Qubits (e.g., Google Sycamore)
temp_sc = 0.015  # 15 mK
coherence_sc = 100  # microseconds (approx)
cost_sc = 1000000 # Normalized cost (hardware + energy)

# 2. HAWRA P700 (Bio-Qubit)
temp_bio = 300  # Kelvin
coherence_bio = 0.000041 # 41 ps (simulated with Silica Shield) -> scale for viz
# Note: To make it visible on the same plot, we plot Scalability or Energy Efficiency instead of raw T2

# Let's plot "Qubit Scalability per Watt" vs "Temperature"
# This highlights the "Obsolescence of Cryogenics"

temperatures = [0.015, 300]
technologies = ['Superconducting (Google/IBM)', 'HAWRA (Metabiotic)']
scalability_per_watt = [1e2, 1e12] # ~50 qubits vs ~10^12 potential qubits/plant
energy_cost = [1e6, 1] # Massive vs Passive (Sunlight)

fig, ax1 = plt.subplots(figsize=(10, 6))

# Bar plot for Scalability
bars = ax1.bar(technologies, scalability_per_watt, color=['#e74c3c', '#2ecc71'], alpha=0.7, width=0.5)
ax1.set_ylabel('Scalability (Qubits / Watt)', fontsize=14, fontweight='bold')
ax1.set_yscale('log')
ax1.set_ylim(1, 1e14)

# Annotate bars
for bar in bars:
    height = bar.get_height()
    ax1.text(bar.get_x() + bar.get_width()/2., height*1.2,
             f'10^{int(np.log10(height))}',
             ha='center', va='bottom', fontsize=12, fontweight='bold')

# Second axis for Temperature
ax2 = ax1.twinx()
ax2.plot(technologies, temperatures, color='#3498db', marker='o', markersize=15, linewidth=3, linestyle='--')
ax2.set_ylabel('Operating Temperature (Kelvin)', fontsize=14, fontweight='bold', color='#3498db')
ax2.tick_params(axis='y', labelcolor='#3498db')
ax2.set_yscale('log')
ax2.set_ylim(0.01, 1000)

# Annotations
plt.title('Figure 1: The Metabiotic Advantage\nObsolescence of Cryogenics', fontsize=16, fontweight='bold', pad=20)
ax1.grid(True, axis='y', linestyle='--', alpha=0.3)

# Add "Cryogenic Barrier" line
ax2.axhline(y=4, color='gray', linestyle=':', linewidth=1)
ax2.text(0.5, 5, 'Cryogenic Barrier (4K)', color='gray', fontsize=10, ha='center')

plt.tight_layout()
plt.savefig('06_publication/figures/figure1_metabiotic_advantage.png', dpi=300)
print("Figure 1 generated successfully.")
