import numpy as np
import matplotlib.pyplot as plt
from lib import detector as d
from lib import particles as p
from lib import experiments as e
from lib import source as s

# File path to save the output spectrum plot
file_path = "/mnt/c/Users/User/Desktop/info/Compton/Simulation/plots/"

def plot_energy_spectrum(energies, fileNamePNG, bins=100, title="Energy Spectrum"):
    """
    Plots a histogram of the detected photon energies.
    
    :param energies: List of detected photon energies in keV.
    :param fileNamePNG: File path to save the plot as a PNG image.
    :param bins: Number of bins for the histogram.
    :param title: Title of the plot.
    """
    # Filter out energies less than or equal to 1 keV, which might be noise or irrelevant
    energies_def = [energy for energy in energies if energy > 1]

    if not energies_def:
        print("No energies to plot.")
        return

    # Create a histogram of detected energies
    plt.figure(figsize=(10, 6))
    plt.hist(energies_def, bins=bins, color="blue", alpha=0.7, edgecolor="blue")
    
    # PP ---> Photopeak, CE ---> Compton Edge
    # Mark the photopeak (PP) and Compton Edge (CE) for the energies of interest
    plt.axvline(511, color='red', linestyle='--', linewidth=1.5, label='PP for 511 keV')
    plt.axvline(1274, color='red', linestyle='--', linewidth=1.5, label='PP for 1274 keV')
    plt.legend(loc="upper right")

    # Set plot title and labels
    plt.title(title, fontsize=16)
    plt.xlabel("Energy (keV)", fontsize=14)
    plt.ylabel("Counts", fontsize=14)
    plt.xlim((0, 900))
    plt.ylim((0, 100))
    plt.grid(alpha=0.4)  # Grid for readability
    plt.tight_layout()  # Adjust layout to fit everything
    plt.savefig(fileNamePNG)  # Save the plot as PNG
    plt.close()  # Close the plot to free up resources

# Initialize detectors with their respective positions and dimensions
franco = d.Detector([0, 45, 0], 2.54, 5.08, 0.1194)  # Detector "Franco"

number_of_photons = 100000  # Number of photons to simulate
source = s.Source({380: 0.8, 511: 0.16, 1275: 0.04})

# Simulate spectroscopy measurements and plot the energy spectrum 
energies = e.spectroscopy_measurement(number_of_photons, franco, source, True)
plot_energy_spectrum(energies, "try.png", 10000)


print("End")

