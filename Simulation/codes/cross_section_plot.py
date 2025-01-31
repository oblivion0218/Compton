import numpy as np
import matplotlib.pyplot as plt
import libraries.particles as p
import libraries.interactions as i
import libraries.target as t

# File path to save the output plots
file_path = "/mnt/c/Users/User/Desktop/info/Compton/Simulation/plots/cross_sections/"

# Define the energy range for plotting (in keV)
energies = np.linspace(10, 1500, 1000)  # Generates 1000 points between 10 keV and 1500 keV

# Define material properties for NaI (sodium iodide)
target = t.Target(
    position=[0, 0, 0],    # Position of the target (arbitrary here)
    radius=5,              # Radius of the target (arbitrary value)
    width=1,               # Width of the target (arbitrary value)
    Z=49.7,                # Effective atomic number for NaI
    density=3.67,          # Density of NaI in g/cm^3
    molar_mass=149.89      # Molar mass of NaI in g/mol
)

# Initialize arrays to store computed values
cross_section_pe = []       # Photoelectric cross-section
cross_section_com = []      # Compton scattering cross-section
attenuation_factors = []    # Attenuation factors for each energy
interaction_prob = []       # Interaction probabilities for each energy

# Distance the photon travels through the material (example: 1 cm)
distance = 1.0

# Compute cross-sections, attenuation factors, and interaction probabilities
for E in energies:
    # Create a photon object with energy E and a default direction
    photon = p.Photon(energy=E, direction=[0, 1, 0])

    # Calculate photoelectric and Compton cross-sections for the photon
    sigma_pe = i.cross_section_photoelectric(photon, target.Z)
    sigma_com = i.cross_section_compton(photon, target.Z)
    
    # Append cross-section values to the respective arrays
    cross_section_pe.append(sigma_pe)
    cross_section_com.append(sigma_com)

    # Calculate interaction probability and attenuation factor
    interaction_prob.append(i.interaction_probability(photon, distance, target))
    attenuation_factors.append(i.attenuation_factor(photon, distance, sigma_pe + sigma_com, target))

# Convert lists to numpy arrays for plotting
cross_section_pe = np.array(cross_section_pe)
cross_section_com = np.array(cross_section_com)
attenuation_factors = np.array(attenuation_factors)
interaction_prob = np.array(interaction_prob)

import matplotlib.pyplot as plt

# Create subplots
fig, axs = plt.subplots(3, 1, figsize=(18, 18), sharex=True)

# Plot cross-sections
axs[0].plot(energies, cross_section_pe, label="Photoelectric Cross Section", color='blue')
axs[0].plot(energies, cross_section_com, label="Compton Cross Section", color='red')
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
axs[2].set_title("Interaction Probability vs Photon Energy")
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

# Compute cross-sections for each Z
for Z in Z_values:
    photon = p.Photon(energy=photon_energy, direction=[0, 1, 0])
    cross_section_pe_vs_Z.append(i.cross_section_photoelectric(photon, Z))
    cross_section_com_vs_Z.append(i.cross_section_compton(photon, Z))

cross_section_pe_vs_Z = np.array(cross_section_pe_vs_Z)
cross_section_com_vs_Z = np.array(cross_section_com_vs_Z)

# Plot cross-sections vs Z
plt.figure(figsize=(12, 8))
plt.plot(Z_values, cross_section_pe_vs_Z, label="Photoelectric Cross Section", color='blue')
plt.plot(Z_values, cross_section_com_vs_Z, label="Compton Cross Section", color='red')

plt.yscale('log')  # Logarithmic scale for cross-sections
plt.xlabel("Atomic Number (Z)")
plt.ylabel("Cross Section " + r"$(cm^2)$")
plt.title("Cross Sections vs Atomic Number (Z) for Energy = 511 keV")
plt.legend()
plt.grid(True, which="both", linestyle="--", linewidth=0.5)
plt.savefig(file_path + "cross_sections_vs_Z.png")  # Save the plot as a PNG file

print("\nEnd\n")