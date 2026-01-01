import numpy as np

t = np.linspace(0, 48, 480)
day = (np.sin(np.pi * t / 24) + 1) / 2
night_atp = 0.35
atp = day * 0.8 + night_atp

data = np.vstack((t, atp)).T
np.savetxt("24h_cam.csv", data, delimiter=",", header="time,atp")

print("Données de production d'ATP sur 24h générées dans 24h_cam.csv")
