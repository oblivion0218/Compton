import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from tqdm import tqdm
from iminuit import Minuit
from iminuit.cost import LeastSquares

N_tot = 1000000  # Total number of photons in each simulation

# file_path = "/mnt/c/Users/User/Desktop/info/Compton/Simulation/simulated_events/"
file_path = "/mnt/c/Users/User/Desktop/info/Gamma-simulation/simulated_events/"

sim_runs = []
photons_left_target = []
interacting_photons = []
photoelectric_interactions = []
compton_1st = []
compton_2nd = []
compton_3rd = []
photons_reached_detector = []
photons_observed = []
energies = []
true_energies = []

file_names = []

if os.path.isdir(file_path + "data/"):
    # Iterate through the items in the directory
    for file_name in os.listdir(file_path + "data/"):
        # Create the full file path
        full_file_path = os.path.join(file_path + "data/", file_name)
        # Add the file name if it is a file (not a directory)
        if os.path.isfile(full_file_path):
            file_names.append(file_path + "data/" + file_name)

# Process each file
for file_name in file_names:
    # Extract run number from filename
    run_number = int(os.path.basename(file_name).split("_")[0])
    sim_runs.append(run_number)
    
    # Parse the metrics from the file
    with open(file_name, 'r') as f:
        lines = f.readlines()
        
        # Extract metrics
        photons_left_target.append(int(lines[0].split(": ")[1]))
        interacting_photons.append(int(lines[1].split(": ")[1]))
        photoelectric_interactions.append(int(lines[2].split(": ")[1]))
        compton_1st.append(int(lines[3].split(": ")[1]))
        compton_2nd.append(int(lines[4].split(": ")[1]))
        compton_3rd.append(int(lines[5].split(": ")[1]))
        photons_reached_detector.append(int(lines[6].split(": ")[1]))
        photons_observed.append(int(lines[7].split(": ")[1]))
        
        capture_true_energies = False
        for i in range(9, len(lines)):
            line = lines[i].strip()
            if line == "True Detected Energies:":
                capture_true_energies = True
                continue
            if line and not capture_true_energies:  # Detected energies
                energies.append(float(line))
            elif line and capture_true_energies:  # True energies
                true_energies.append(float(line))

# Convert lists to numpy arrays for easier manipulation
sim_runs = np.array(sim_runs)
photons_left_target = np.array(photons_left_target)
interacting_photons = np.array(interacting_photons)
photoelectric_interactions = np.array(photoelectric_interactions)
compton_1st = np.array(compton_1st)
compton_2nd = np.array(compton_2nd)
compton_3rd = np.array(compton_3rd)
photons_reached_detector = np.array(photons_reached_detector)
photons_observed = np.array(photons_observed)
energies = np.array(energies)
true_energies = np.array(true_energies)

plt.figure(figsize=(12, 8))

# Plotting all metrics on the same plot
plt.plot(sim_runs, interacting_photons/N_tot, 'o', color='red', label='Photon Interactions')
plt.plot(sim_runs, photoelectric_interactions/N_tot, 'o', color='green', label='Photoelectric Interactions')
plt.plot(sim_runs, compton_1st/N_tot, 'o', color='purple', label='Compton 1st')
plt.plot(sim_runs, compton_2nd/N_tot, 'o', color='orange', label='Compton 2nd')
plt.plot(sim_runs, compton_3rd/N_tot, 'o', color='brown', label='Compton 3rd')

plt.title('Photon Interaction probability across simulation runs')
plt.xlabel('Simulation Run')
plt.ylabel('Probability / Scaled Efficiency')
plt.legend(fontsize=12)
plt.grid(True)
plt.savefig(file_path + "plots/interaction_probabilities.png")

plt.figure(figsize=(12, 8))
plt.plot(sim_runs, photons_reached_detector/N_tot, 'o', color='red', label='Photons Reaching Detector')
plt.plot(sim_runs, photons_observed/N_tot, 'o', color='blue', label='Photons Observed')
plt.title('Photon interaction across first simulation runs')
plt.xlabel('Simulation run')
plt.ylabel('Probability / Scaled Efficiency')
plt.legend(fontsize=12)
plt.grid(True)
plt.savefig(file_path + "plots/detector_parameter.png")

detection_efficiency = photons_observed / photons_reached_detector
sorted_indices = np.argsort(sim_runs)
sim_runs_sorted = sim_runs[sorted_indices]
detection_efficiency_sorted = detection_efficiency[sorted_indices]

plt.figure(figsize=(12, 8))
plt.plot(sim_runs, detection_efficiency / 100, 'o', color='black', label='Detection Efficiency (scaled)')
plt.title('Detection Efficiency across simulation runs')
plt.xlabel('Simulation Run')
plt.ylabel('Detection Efficiency (scaled)')
plt.grid(True)
plt.savefig(file_path + "plots/detection_efficiency.png")

# Create a new combined histogram with Freedman-Diaconis binning
plt.figure(figsize=(12, 8))

# Create histograms
plt.hist(energies, bins=100, alpha=0.4, color='blue', 
         edgecolor='black', label='Detected Energies')
plt.hist(true_energies, bins=100, color='red', 
              histtype="step", label='True Energies')

plt.title('Probability Density of Photon Energies')
plt.xlabel('Energy (keV)')
plt.ylabel('Counts')
plt.legend(fontsize=12)
plt.grid(True)
plt.savefig(file_path + "plots/photon_energies.png")


