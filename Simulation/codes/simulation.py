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
    plt.axvline(341, color='orange', linestyle='--', linewidth=1.5, label='CE for 511 keV') # 2/3 511
    plt.axvline(1061, color='orange', linestyle='--', linewidth=1.5, label='CE for 1274 keV') 
    plt.legend(loc="upper right")

    # Set plot title and labels
    plt.title(title, fontsize=16)
    plt.xlabel("Energy (keV)", fontsize=14)
    plt.ylabel("Counts", fontsize=14)
    plt.yscale("log")  # Log scale for better visibility of data
    plt.grid(alpha=0.4)  # Grid for readability
    plt.tight_layout()  # Adjust layout to fit everything
    plt.savefig(fileNamePNG)  # Save the plot as PNG
    plt.close()  # Close the plot to free up resources

def visualization_3D(fileNamePNG, detectors, photons, target=None):
    """
    Creates a 3D visualization of the photon hit points and detector positions.
    
    :param fileNamePNG: Path to save the output 3D visualization image.
    :param detectors: List of detector objects.
    :param photons: List of photon objects, each with a position attribute.
    :param target: Optional target object to visualize.
    """
    hit_points = np.array([photon.position for photon in photons])  # Extract photon hit positions
    
    # Create figure with a single 3D plot
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    # Plot hit points
    ax.scatter(hit_points[:, 0], hit_points[:, 1], hit_points[:, 2], 
               s=10, color='red', marker='x', label='Hit Points')
    
    # Draw detectors and target
    [detector.draw_3D(ax) for detector in detectors]
    if target:
        target.draw_3D(ax)
    
    # Add labels and title
    ax.set_xlabel('X (cm)', fontsize=12)
    ax.set_ylabel('Y (cm)', fontsize=12)
    ax.set_zlabel('Z (cm)', fontsize=12)
    ax.set_title('3D Photon Hit Positions', fontsize=14)
    
    # Improve visibility
    ax.legend(loc='upper right', fontsize=10)
    ax.grid(True, alpha=0.3)
    
    # Set equal aspect ratio for better visualization
    ax.set_box_aspect([1, 1, 1])
    
    # Adjust layout and save
    plt.tight_layout()
    plt.savefig(fileNamePNG, dpi=300)
    plt.close()

# Initialize detectors with their respective positions and dimensions
ugo = d.Detector(([0, 15, 0], [0, 17.58, 0]), 1.27, 0.0903)  # Detector "Ugo"
franco = d.Detector(([0, -15, 0], [0, -20.08, 0]), 2.54, 0.0695)  # Detector "Franco"

number_of_photons = 1000000  # Number of photons to simulate

# Simulate spectroscopy measurements and plot the energy spectrum 
# energies = e.spectroscopy_measurement(number_of_photons, ugo, s.Source(), True)
# plot_energy_spectrum(energies, file_path + "spectrum/spettroscopy_ugo.png", 10000)
# energies = e.spectroscopy_measurement(number_of_photons, franco, s.Source(), True)
# plot_energy_spectrum(energies, file_path + "spectrum/spettroscopy_franco.png", 10000)
# print("Spectroscopy measurements completed.")

# Simulate photon emission from a source
source = s.Source()  # Create a source object
photons = source.photon_emission(1000, np.pi/4, np.pi/6)  # Simulate photon emission
[photon.propagation(np.linalg.norm(franco.position)) for photon in photons]  # Propagate the photons

franco.rotate(np.pi / 4, [0, 0, 0], "z")
# Visualize the 3D photon hit positions and detectors for coincidence measurement
visualization_3D(file_path + "3D_plots/3D_visualization.png", [ugo, franco], photons)
print("3D visualization completed.")

print("End")

