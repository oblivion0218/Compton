import numpy as np
import matplotlib.pyplot as plt
import random
import particles as p
import detector as d
import source as s
import interactions as i 

file_path = "/mnt/c/Users/User/Desktop/info/Compton/Simulation/"

ugo = d.Detector([0, 5, 0], 1.27, 2.58, 0.0695)
franco = d.Detector([0, -5, 0], 2.54, 5.08, 0.903)
source = s.Source()

def gamma_detection(photons: list[p.Photon], detector: d.Detector, step: float = 0.1, source: s.Source = source):
    distance = np.linalg.norm(detector.position - source.position)

    detected_energies = []

    for photon in photons: 
        photon.propagation(distance)

        alpha_angle = np.arctan(photon.position[0] / photon.position[2]) if (photon.position[0] != 0 and photon.position[2] != 0) else 0
        h = distance - photon.position[1]  
        l = h / np.cos(alpha_angle)

        distance += l
        photon.propagation(distance)

        r = 0
        while detector.is_in_detector(photon.position):
            if random.uniform(0, 1) < i.interaction_probability(photon, r):
                interaction = i.Interaction(i.which_interaction(photon, detector))
                electron = interaction.interaction(photon)

                detected_energies.append(detector.detection(electron))
                break

            r += step
            photon.propagation(distance + step)
    
    return detected_energies

def plot_energy_spectrum(energies, fileNamePNG, bins=100, title="Energy Spectrum"):
    """Plot a histogram of the energy spectrum."""
    if not energies:
        print("No energies to plot.")
        return

    plt.figure(figsize=(10, 6))
    plt.hist(energies, bins=bins, color="blue", alpha=0.7, edgecolor="blue")
    
    plt.axvline(511, color='red', linestyle='--', linewidth=1.5, label='511 keV')
    plt.axvline(1274, color='red', linestyle='--', linewidth=1.5, label='1274 keV')
    plt.legend()

    plt.title(title, fontsize=16)
    plt.xlabel("Energy (keV)", fontsize=14)
    plt.ylabel("Counts", fontsize=14)
    plt.ylim(0, 2500)
    plt.grid(alpha=0.4)
    plt.tight_layout()
    plt.savefig(fileNamePNG)
    plt.close()

number_of_photons = 100000
energies = gamma_detection(source.testing_photons(number_of_photons), ugo)
plot_energy_spectrum(energies, file_path + "Simulation.png", 200)
