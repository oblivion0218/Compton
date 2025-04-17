import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from lib import particles as p
from scipy.integrate import quad
import random

# File path to save the output spectrum plot
# file_path = "/home/leonardo/Compton/Simulation/plots/compton_angles_distributions/" # LEO
file_path = "/mnt/c/Users/User/Desktop/info/Compton/Simulation/plots/compton_angles_distributions/"  # ANDRE

# Instantiate a Photon object with arbitrary energy and direction

number_of_photons = 100000
photon = p.Photon(511, [0, 0, 1])

"""
Klein-Nishina distibution
"""

# Generate 10,000 random Compton scattering angles
random_angles = [photon.compton_angle() for _ in range(number_of_photons)]

# Define the Klein-Nishina distribution over a range of angles
angle_range = np.linspace(0, np.pi, 500)

# Compute the unnormalized Klein-Nishina PDF values
pdf_values = [photon.klein_nishina(angle) for angle in angle_range]

# Integrate the Klein-Nishina function over the full range to normalize
normalization_factor, _ = quad(photon.klein_nishina, 0, np.pi)

# Normalize the PDF values
normalized_pdf_values = [value / normalization_factor for value in pdf_values]

# Generate 10,000 random electron energies
random_electron_energies = [photon.energy - photon.compton_scattering(angle) for angle in random_angles]

# Plot histogram of random angles, normalized Klein-Nishina PDF, and electron energy distribution in the same figure
plt.figure(figsize=(20, 8))

plt.subplot(1, 2, 1)
# Histogram of random angles
plt.hist(random_angles, bins=500, density=True, alpha=0.7, label='Random Angles')
# Normalized Klein-Nishina PDF
plt.plot(angle_range, normalized_pdf_values, color='red', label='Normalized Klein-Nishina PDF')

plt.xlabel('Scattering Angle (radians)')
plt.ylabel('Probability Density')
plt.title('Compton Scattering: Angle Distributions')
plt.legend()
plt.grid(True)

plt.subplot(1, 2, 2)
# Histogram of random electron energies
plt.hist(random_electron_energies, bins=1000, alpha=0.7, label='Random Electron Energies', color='green')

plt.xlabel('Electron Energy (keV)')
plt.ylabel('Counts')
plt.title('Compton Scattering: Electron Energy Distributions (511 keV Photons)')
plt.legend()
plt.grid(True)

# Save and show the figure
plt.savefig(file_path + 'angle_distributions_Klein_Nishina.png')


norm = mcolors.Normalize(vmin=np.min(random_electron_energies), vmax=np.max(random_electron_energies))
cmap = cm.viridis

plt.figure(figsize=(10, 8))
ax = plt.subplot(1, 1, 1, polar=True)
photon_energies_range = [photon.energy - photon.compton_scattering(angle) for angle in angle_range]


# Disegna la curva segmentando e colorando ogni tratto in base al valore di energia
for i in range(len(angle_range) - 1):
    theta_seg = angle_range[i:i+2]
    r_seg = normalized_pdf_values[i:i+2]
    # Il colore del segmento viene scelto in base all'energia del punto corrente
    color = cmap(norm(photon_energies_range[i]))
    ax.plot(theta_seg, r_seg, color=color, linewidth=2)

# Aggiungi la colorbar per indicare le energie
sm = cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax, pad=0.1)
cbar.set_label('Energy of the scattered electron (keV)')

# Optional: prettify
ax.set_theta_zero_location("N")  # 0 rad at the top
ax.set_theta_direction(-1)       # clockwise
ax.set_title('Klein-Nishina Distribution (Radial Plot)', va='bottom')
ax.grid(True)

# Save the figure
plt.savefig(file_path + 'angle_distributions_Klein_Nishina_radial.png')
