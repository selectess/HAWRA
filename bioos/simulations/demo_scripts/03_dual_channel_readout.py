import numpy as np

trials = 10000
p_atp = 0.62
p_ros = 1 - p_atp

green = np.random.choice([0, 1], trials, p=[1 - p_atp, p_atp])
red = np.random.choice([0, 1], trials, p=[1 - p_ros, p_ros])

data = np.vstack((green, red)).T
np.savetxt("dual_channel.csv", data, delimiter=",", header="green,red", fmt="%d")

print("Données de lecture double canal générées dans dual_channel.csv")
