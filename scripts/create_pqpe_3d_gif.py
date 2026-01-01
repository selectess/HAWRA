import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
from Bio import SeqIO

base = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
gb_path = os.path.join(base, '01_genomics', 'plasmids', 'HAWRA_FINAL_VALIDATED.gb')
out_dir = os.path.join(base, '05_data', 'results')
os.makedirs(out_dir, exist_ok=True)
out_gif = os.path.join(out_dir, 'pqpe_plasmid_3d.gif')

record = SeqIO.read(gb_path, 'genbank')
seq_len = len(record.seq)
features = []
for f in record.features:
    try:
        start = int(f.location.start)
        end = int(f.location.end)
    except Exception:
        continue
    features.append({
        'start': start,
        'end': end,
        'type': f.type,
        'label': f.qualifiers.get('gene', [''])[0]
    })

R = 5.0
r = 1.0
theta = np.linspace(0, 2 * np.pi, 1000)
x_center = (R + r * np.cos(0)) * np.cos(theta)
y_center = (R + r * np.cos(0)) * np.sin(theta)
z_center = r * np.sin(0) * np.ones_like(theta)

def arc_indices(start, end):
    a0 = (start / seq_len) * 2 * np.pi
    a1 = (end / seq_len) * 2 * np.pi
    if a1 < a0:
        a0, a1 = a1, a0
    mask = (theta >= a0) & (theta <= a1)
    return np.where(mask)[0]

def feature_color(label, ftype):
    palette = {
        'psaA': '#4b0082',
        'CRY2': '#0000ff',
        'Luc': '#00ff00',
        'Lsi1': '#ffA500',
        'HSP70': '#ff4500',
        'PEPC': '#dc143c',
        'P700': '#ff1493'
    }
    if label in palette:
        return palette[label]
    if ftype == 'CDS':
        return 'cyan'
    return 'grey'

fig = plt.figure(figsize=(10, 5))
ax3d = fig.add_subplot(1, 2, 1, projection='3d')
ax2 = fig.add_subplot(1, 2, 2)
ax2.axis('off')
ax2.text(0.0, 0.9, 'HAWRA', fontsize=14, weight='bold')
ax2.text(0.0, 0.8, 'Simulateur unifié: Environnement + Biologie + Quantique', fontsize=10)
ax2.text(0.0, 0.65, 'Arbol → BSIM', fontsize=14, weight='bold')
ax2.text(0.0, 0.55, 'Langage expérimental → Bytecode JSON', fontsize=10)
ax2.text(0.0, 0.4, 'BioOS', fontsize=14, weight='bold')
ax2.text(0.0, 0.3, 'Orchestrateur: compile, exécute, enregistre', fontsize=10)
ax2.text(0.0, 0.15, 'Programmabilités', fontsize=14, weight='bold')
ax2.text(0.0, 0.05, 'Stimuli, GRN, circuits quantiques, séquencement', fontsize=10)

ax3d.plot(x_center, y_center, z_center, color='black', linewidth=1.5)
for f in features:
    idx = arc_indices(f['start'], f['end'])
    if idx.size == 0:
        continue
    ax3d.plot(x_center[idx], y_center[idx], z_center[idx], color=feature_color(f['label'], f['type']), linewidth=3)

ax3d.set_xlim(-R - r - 1, R + r + 1)
ax3d.set_ylim(-R - r - 1, R + r + 1)
ax3d.set_zlim(-r - 1, r + 1)
ax3d.set_xticks([])
ax3d.set_yticks([])
ax3d.set_zticks([])
ax3d.set_title('PQPE Plasmide (3D)')

def animate(i):
    angle = i * (360.0 / 60)
    ax3d.view_init(elev=20, azim=angle)
    return fig

anim = FuncAnimation(fig, animate, frames=60, interval=50, blit=False)
try:
    anim.save(out_gif, writer='imagemagick', fps=12)
except Exception:
    from matplotlib.animation import PillowWriter
    anim.save(out_gif, writer=PillowWriter(fps=12))
print(f'Sauvegardé: {out_gif}')