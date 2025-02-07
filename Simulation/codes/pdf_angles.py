import numpy as np
import matplotlib.pyplot as plt
from lib import particles as p
from scipy.integrate import quad
import random

# File path to save the output spectrum plot
file_path = "/mnt/c/Users/User/Desktop/info/Compton/Simulation/plots/compton_angles_distributions/"

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

"""
Uniform distribution
"""

# Generate 10,000 random Compton scattering angles (uniform distribution)
random_angles = [random.uniform(0, np.pi) for _ in range(number_of_photons)]

# Define the uniform PDF over a range of angles
angle_range = np.linspace(0, np.pi, 500)
uniform_pdf_values = [1 / np.pi for _ in angle_range]  # Uniform distribution over [0, pi]

# Generate 10,000 random electron energies
random_electron_energies = [photon.energy - photon.compton_scattering(angle) for angle in random_angles]

# Plot histogram of random angles, uniform PDF, and electron energy distribution in the same figure
plt.figure(figsize=(20, 8))

plt.subplot(1, 2, 1)
# Histogram of random angles
plt.hist(random_angles, bins=500, density=True, alpha=0.7, label='Random Angles')
# Uniform PDF
plt.plot(angle_range, uniform_pdf_values, color='red', label='Uniform PDF')

plt.xlabel('Scattering Angle (radians)')
plt.ylabel('Probability Density')
plt.title('Compton Scattering: Angle Distributions (Uniform PDF)')
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
plt.savefig(file_path + 'angle_distributions_uniform.png')