import geometric_simulation as gs
import numpy as np
import matplotlib.pyplot as plt

path_file = "/mnt/c/Users/User/Desktop/info/Compton/Simulation/"

activity = 130000  # Photons per second
n_photons = 10000  # Number of photons simulated for visualization

# Generate random directions
directions = gs.random_directions(n_photons)

# # Visualize hits, projections, and detection counts for Ugo and Franco detectors
gs.visualize_coincidence(path_file + "coincidence.png", directions)

n_tries = 1000
n_photons = 1000000
distance = 5
ugo_hits = []
franco_hits = []

for i in range(n_tries):
    directions = gs.random_directions(n_photons)
    hit_points = np.array([gs.hit_position(dir, distance) for dir in directions])
    ugo_hits.append(gs.is_in_detector(hit_points, "ugo"))
    franco_hits.append(gs.is_in_detector(hit_points, "franco"))

plt.figure(figsize=(12, 6))

# Histogram for Ugo hits
plt.subplot(1, 2, 1)
plt.hist(ugo_hits, bins=200, color='blue', alpha=0.7, edgecolor='blue')
plt.title("Histogram of Hits in Ugo Detector")
plt.xlabel("Number of Hits")
plt.ylabel("Frequency")

# Histogram for Franco hits
plt.subplot(1, 2, 2)
plt.hist(franco_hits, bins=200, color='orange', alpha=0.7, edgecolor='orange')
plt.title("Histogram of Hits in Franco Detector")
plt.xlabel("Number of Hits")
plt.ylabel("Frequency")

# Show the histograms
plt.tight_layout()
plt.savefig(path_file + "histograms.png")