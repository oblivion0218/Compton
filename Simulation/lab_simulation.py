import geometric_simulation as gs

path_file = "/mnt/c/Users/User/Desktop/info/Compton/Simulation/"

activity = 130000  # Photons per second
n_photons = 5000  # Number of photons simulated for visualization

# Generate random directions
directions = gs.random_directions(n_photons)

# # Visualize hits, projections, and detection counts for Ugo and Franco detectors
gs.visualize_coincidence(path_file + "coincidence.png", directions)