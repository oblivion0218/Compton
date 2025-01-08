import numpy as np
import matplotlib.pyplot as plt
import particles as p
import detector as d
import interactions as i

# File path to save the output spectrum plot
file_path = "/mnt/c/Users/User/Desktop/info/Gamma-simulation/plots/"

# Energy range for plotting (in keV)
energies = np.linspace(10, 1500, 1000)

detector = d.Detector([0, 15, 0], 1.27, 2.58, 0.0903, Z=49.7)

# Arrays to store cross sections
cross_section_pe = []
cross_section_com = []

# Compute cross sections for each energy
for E in energies:
    photon = p.Photon(energy=E, direction=[0, 1, 0])
    sigma_pe = i.cross_section_photoelectric(photon, detector)
    sigma_com = i.cross_section_compton(photon, detector)
    
    cross_section_pe.append(sigma_pe)
    cross_section_com.append(sigma_com)

# Convert to numpy arrays for plotting
cross_section_pe = np.array(cross_section_pe)
cross_section_com = np.array(cross_section_com)

# Plotting
plt.figure(figsize=(10, 6))
plt.plot(energies, cross_section_pe, label="Photoelectric Cross Section", color='blue')
plt.plot(energies, cross_section_com, label="Compton Cross Section", color='red')
plt.xscale('log')
plt.yscale('log')
plt.xlabel("Photon Energy (keV)")
plt.ylabel("Cross Section (m²)")
plt.title("Cross Sections vs Photon Energy")
plt.legend()
plt.grid(True, which="both", linestyle="--", linewidth=0.5)
plt.savefig(file_path + "cross_sections.png")
