import numpy as np
import matplotlib.pyplot as plt
import random
from tqdm import tqdm
from lib import interactions as i
from lib import detector as d
from lib import experiments as e
from lib import source as s
from lib import visualization as v

# Object initialization
# *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
source = s.Source({511: 1, 1274: 0})  # Create a source object
target = d.Target(([0, 5, 0], [0, 6, 0]), 3)  # Create a target object
detector = d.Detector(([0, 33.04, 0], [0, 35.58, 0]), 2.54, 0.0695)  # Detector "Franco"

# Initial parameters
# *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
step = 0.025
number_of_photons = 1000000

theta = 35

r_gate = 1.27
d_gate_source = 18.54
theta_gate = np.arctan(r_gate/d_gate_source)

# Simulation of photons emitted from the source 
# *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
photons = source.photon_emission(number_of_photons, theta_gate, 2 * np.pi, axis="y", forward_backward=False) 
theta_rad = np.radians(theta)  # Convert theta to radians
# target_theta = theta_rad + (np.pi - theta_rad) / 2
detector.rotate(theta_rad, [0, 5.5, 0], "z")

print(detector.position[0])

# target.rotate(target_theta, [0, 5.5, 0], "z")

[e.photon_propagation_to_target(photon, target) for photon in photons] 
# v.visualization_3D_plotly("photons_" + str(theta) + ".html", [detector], photons, source=source, target=target)

ip = []

# Photons interacting with the target
# *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
for photon in tqdm(photons, desc="Simulating photons", unit="photon"):
    photon.propagation(step)
    while target.is_inside(photon.position):
        if i.interaction_probability(photon, step, target) > random.uniform(0, 1):
            ip.append(photon.position)
            break
        else:
            photon.propagation(step)

# Save the results in a txt file
with open("IP_simulation_" + str(theta) + ".txt", "a") as f:
    for l in ip:
        f.write(str(l) + "\n")
f.close()

