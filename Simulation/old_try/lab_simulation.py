import geometric_simulation as gs
import numpy as np
import matplotlib.pyplot as plt

path_file = "/mnt/c/Users/User/Desktop/info/Compton/Simulation/"

activity = 130000  # Photons per second
n_photons = 10000  # Number of photons simulated for visualization

directions = gs.random_directions(n_photons)
gs.visualize_coincidence(path_file + "coincidence_5cm.png", directions, 5)

# Simulation parameters
n_tries = 1000
distances = [5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
mean_ugo_hits = []
mean_franco_hits = []

# Loop through distances
for distance in distances:
    ugo_hits = []
    franco_hits = []

    # Perform trials for each distance
    for _ in range(n_tries):
        directions = gs.random_directions(n_photons)
        hit_points = np.array([gs.hit_position(dir, distance) for dir in directions])
        ugo_hits.append(np.sum(gs.is_in_detector(hit_points, "ugo")) / n_photons)
        franco_hits.append(np.sum(gs.is_in_detector(hit_points, "franco")) / n_photons)

    # Calculate means for the current distance
    mean_ugo_hits.append(np.mean(ugo_hits))
    mean_franco_hits.append(np.mean(franco_hits))

    # Plot histograms for the current distance
    plt.figure(figsize=(14, 6))


    # Histogram for Ugo hits
    plt.subplot(1, 2, 1)
    ugo_mean = np.mean(ugo_hits) 
    plt.hist(ugo_hits, bins=100, color='blue', alpha=0.7, edgecolor='blue')
    plt.title(f"Histogram of Hits in Ugo Detector (Distance: {distance} cm)")
    plt.xlabel("Number of Hits / Number of Photons")
    plt.ylabel("Frequency")

    plt.text(plt.gca().get_xlim()[1] * 0.95, plt.gca().get_ylim()[1] * 0.95, # Posizione del testo
             f"Mean: {ugo_mean:.4f}", fontsize=10, color='blue', bbox=dict(facecolor='white', alpha=0.8))

    # Histogram for Franco hits
    plt.subplot(1, 2, 2)
    franco_mean = np.mean(franco_hits)  
    plt.hist(franco_hits, bins=100, color='orange', alpha=0.7, edgecolor='orange')
    plt.title(f"Histogram of Hits in Franco Detector (Distance: {distance} cm)")
    plt.xlabel("Number of Hits / Number of Photons")
    plt.ylabel("Frequency")

    plt.text(plt.gca().get_xlim()[1] * 0.95, plt.gca().get_ylim()[1] * 0.95, # Posizione del testo
             f"Mean: {franco_mean:.4f}", fontsize=10, color='orange', bbox=dict(facecolor='white', alpha=0.8))

    # Save histograms
    plt.tight_layout()
    plt.savefig(path_file + f"Histograms/histograms_{distance}cm_{n_photons}photons.png")
    plt.close()

# Plot mean hits vs distances
plt.figure(figsize=(10, 6))
plt.plot(distances, mean_ugo_hits, label="Ugo Detector", marker='o', color='blue')
plt.plot(distances, mean_franco_hits, label="Franco Detector", marker='s', color='orange')
plt.title(f"Mean Hits vs. Distance for {n_photons} Photons")
plt.xlabel("Distance (cm)")
plt.ylabel(f"Mean Number of Hits of {n_tries} Trials")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig(path_file + "mean_hits_vs_distances.png")
