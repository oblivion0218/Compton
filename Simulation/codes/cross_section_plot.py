import numpy as np
import matplotlib.pyplot as plt
from lib import particles as p
from lib import interactions as i
from lib import detector as d
from tqdm import tqdm

# File path to save the output plots
file_path = "/mnt/c/Users/User/Desktop/info/Gamma-simulation/plots/Interaction_probability/"

# Define the energy range for plotting (in keV)
energies = np.linspace(10, 1500, 1000)  # Generates 1000 points between 10 keV and 1500 keV

# Define material properties for NaI, 1x1 inch target
target = d.Target(
    position=([0, 0, 0], [0, 2.54, 0]),    # Position of the target (arbitrary here)
    radius=1,              # Radius of the target (arbitrary value)
    Z=49.7,                # Effective atomic number for NaI
    density=3.67,          # Density of NaI in g/cm^3
    molar_mass=149.89      # Molar mass of NaI in g/mol
)

# Initialize arrays to store computed values
cross_section_pe = []       # Photoelectric cross-section
cross_section_com = []      # Compton scattering cross-section
total_cross_section = []      # Total cross-section
attenuation_factors = []    # Attenuation factors for each energy
interaction_prob = []       # Interaction probabilities for each energy

# Compute cross-sections, attenuation factors, and interaction probabilities
for E in tqdm(energies, desc="Energies", unit="keV"):
    # Create a photon object with energy E and a default direction
    photon = p.Photon(energy=E, direction=[0, 1, 0])

    # Calculate photoelectric and Compton cross-sections for the photon
    sigma_pe = i.cross_section_photoelectric(photon, target.Z)
    sigma_com = i.cross_section_compton(photon, target.Z)
    sigma_tot = sigma_pe + sigma_com
    
    # Append cross-section values to the respective arrays
    cross_section_pe.append(sigma_pe)
    cross_section_com.append(sigma_com)
    total_cross_section.append(sigma_tot)

    attenuation_factor = i.attenuation_factor(sigma_tot, target)

    # Calculate interaction probability and attenuation factor
    interaction_prob.append(i.interaction_probability(photon, 0.1, target))  
    attenuation_factors.append(attenuation_factor)

# Convert lists to numpy arrays for plotting
cross_section_pe = np.array(cross_section_pe)
cross_section_com = np.array(cross_section_com)
total_cross_section = np.array(total_cross_section)
attenuation_factors = np.array(attenuation_factors)
interaction_prob = np.array(interaction_prob)

import matplotlib.pyplot as plt

# Create subplots
fig, axs = plt.subplots(3, 1, figsize=(18, 18), sharex=True)

# Upper title
fig.suptitle("Cross Sections and Interaction Probability vs Photon Energy for NaI detector with width of 0.1 cm\n", fontsize=12)
# Plot cross-sections
axs[0].plot(energies, cross_section_pe, label="Photoelectric Cross Section", color='blue')
axs[0].plot(energies, cross_section_com, label="Compton Cross Section", color='red')
axs[0].plot(energies, total_cross_section, label="Total Cross Section", color='black', linestyle='--')
axs[0].set_yscale('log')  # Logarithmic scale for cross-section
axs[0].set_ylabel("Cross Section " + r"$(cm^2)$")
axs[0].set_title("Cross Sections vs Photon Energy")
axs[0].legend()
axs[0].grid(True, which="both", linestyle="--", linewidth=0.5)  # Add grid for better visualization

# Plot attenuation factor
axs[1].plot(energies, attenuation_factors, color='green')
axs[1].set_yscale('log')  # Logarithmic scale for attenuation factor
axs[1].set_ylabel("Attenuation Factor " + r"$(cm^{-1})$")
axs[1].set_title("Attenuation Factor vs Photon Energy")
axs[1].grid(True, which="both", linestyle="--", linewidth=0.5)

# Plot interaction probability
axs[2].plot(energies, interaction_prob, color='purple')
axs[2].set_yscale('log')  # Logarithmic scale for interaction probability
axs[2].set_xlabel("Photon Energy (keV)")
axs[2].set_ylabel("Interaction Probability")
axs[2].set_title("Interaction Probability vs Photon Energy for distance = 0.1 cm")
axs[2].grid(True, which="both", linestyle="--", linewidth=0.5)

# Adjust layout and save the combined figure
plt.tight_layout()
plt.savefig(file_path + "cross_section_&_interaction_probability.png")  # Save the image as a PNG file

# Define a range of Z values
Z_values = np.linspace(1, 100, 100)  # Z from 1 (hydrogen) to 100 (fermium)

# Energy of the photon (e.g., 100 keV)
photon_energy = 511  # keV
cross_section_pe_vs_Z = []
cross_section_com_vs_Z = []
total_cross_section_vs_Z = []

# Compute cross-sections for each Z
for Z in tqdm(Z_values, desc="Z values", unit="Z"):
    photon = p.Photon(energy=photon_energy, direction=[0, 1, 0])
    cross_section_pe_vs_Z.append(i.cross_section_photoelectric(photon, Z))
    cross_section_com_vs_Z.append(i.cross_section_compton(photon, Z))
    total_cross_section_vs_Z.append(cross_section_pe_vs_Z[-1] + cross_section_com_vs_Z[-1])

cross_section_pe_vs_Z = np.array(cross_section_pe_vs_Z)
cross_section_com_vs_Z = np.array(cross_section_com_vs_Z)
total_cross_section_vs_Z = np.array(total_cross_section_vs_Z)

# Plot cross-sections vs Z
plt.figure(figsize=(12, 8))
plt.plot(Z_values, cross_section_pe_vs_Z, label="Photoelectric Cross Section", color='blue')
plt.plot(Z_values, cross_section_com_vs_Z, label="Compton Cross Section", color='red')
plt.plot(Z_values, total_cross_section_vs_Z, label="Total Cross Section", color='black', linestyle='--')

plt.yscale('log')  # Logarithmic scale for cross-sections
plt.xlabel("Atomic Number (Z)")
plt.ylabel("Cross Section " + r"$(cm^2)$")
plt.title("Cross Sections vs Atomic Number (Z) for Energy = 511 keV")
plt.legend()
plt.grid(True, which="both", linestyle="--", linewidth=0.5)
plt.savefig(file_path + "cross_sections_vs_Z.png")  # Save the plot as a PNG file

print("\nEnd\n")