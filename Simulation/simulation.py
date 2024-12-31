import numpy as np
import matplotlib.pyplot as plt
import random
import particles as p
import detector as d
import source as s
import interactions as i 

# File path to save the output spectrum plot
file_path = "/mnt/c/Users/User/Desktop/info/Compton/Simulation/"

# Initialize detectors with their respective positions and dimensions
ugo = d.Detector([0, 5, 0], 1.27, 2.58, 0.0903 )  # Detector "Ugo" at (0, 5, 0)
franco = d.Detector([0, -5, 0], 2.54, 5.08, 0.0695)  # Detector "Franco" at (0, -5, 0)
source = s.Source()

# Function to simulate the detection of gamma photons
def gamma_detection(photons: list[p.Photon], detector: d.Detector, step: float = 0.1, source: s.Source = source):
    """
    Simulates the interaction of photons with a detector and calculates the detected energies.
    
    :param photons: List of Photon objects emitted by the source.
    :param detector: Detector object where photons are being detected.
    :param step: Step size for photon propagation in cm.
    :param source: Source object representing the gamma-ray source.
    :return: List of detected photon energies in keV.
    """
    # Calculate the straight-line distance between the source and the detector
    distance = np.linalg.norm(detector.position - source.position)

    # List to store detected photon energies
    detected_energies = []

    for photon in photons: 
        # Propagate the photon from the source to the detector's surface
        photon.propagation(distance)

        # Calculate the angle (α) between the photon's trajectory and the detector's surface
        alpha_angle = np.arctan(photon.position[0] / photon.position[2]) if (photon.position[0] != 0 and photon.position[2] != 0) else 0
        # Calculate the vertical offset (H) and the diagonal length (L) to the detector
        h = distance - photon.position[1]  
        l = h / np.cos(alpha_angle)

        # Adjust the photon's propagation distance to account for the offset
        distance += l
        photon.propagation(distance)

        # Initialize a distance step for interaction testing
        r = 0
        while detector.is_in_detector(photon.position):  # Check if the photon is still inside the detector

            # while photon.energy > 0:
            # Determine if the photon interacts with the detector
            if random.uniform(0, 1) < i.interaction_probability(photon, r):
                # Simulate the interaction and create an electron based on interaction type
                interaction = i.Interaction(i.which_interaction(photon, detector))
                electron = interaction.interaction(photon)

                # Detect the energy deposited by the electron in the detector
                detected_energies.append(detector.detection(electron))
                break  # Stop further propagation after interaction

            # Increment the step distance and propagate the photon further
            r += step
            photon.propagation(distance + step)
    
    return detected_energies

# Function to plot the energy spectrum of detected photons
def plot_energy_spectrum(energies, fileNamePNG, bins=100, title="Energy Spectrum"):
    """
    Plots a histogram of the detected photon energies.
    
    :param energies: List of detected photon energies in keV.
    :param fileNamePNG: File path to save the plot as a PNG image.
    :param bins: Number of bins for the histogram.
    :param title: Title of the plot.
    """
    if not energies:
        print("No energies to plot.")
        return

    # Create a histogram of detected energies
    plt.figure(figsize=(10, 6))
    plt.hist(energies, bins=bins, color="blue", alpha=0.7, edgecolor="blue")
    
    # Mark the photopeak positions (511 keV and 1274 keV) with vertical red dashed lines
    plt.axvline(511, color='red', linestyle='--', linewidth=1.5, label='511 keV')
    plt.axvline(1274, color='red', linestyle='--', linewidth=1.5, label='1274 keV')
    plt.axvline(340, color='orange', linestyle='--', linewidth=1.5, label='CE for 511 keV') # 2/3 511
    plt.legend()

    plt.title(title, fontsize=16)
    plt.xlabel("Energy (keV)", fontsize=14)
    plt.ylabel("Counts", fontsize=14)
    plt.yscale("log")
    # plt.ylim(0, 2500) 
    plt.grid(alpha=0.4)  
    plt.tight_layout()
    plt.savefig(fileNamePNG)  
    plt.close()

number_of_photons = 100000
energies = gamma_detection(source.testing_photons(number_of_photons), ugo)
plot_energy_spectrum(energies, file_path + "Simulation.png", 200)

print("\nEnd\n")
