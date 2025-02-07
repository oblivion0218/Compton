import numpy as np
import random  
from lib import target as t  
from lib import experiments as e 
from lib import particles as p
from lib import source as s  
from lib import interactions as i 
import matplotlib.pyplot as plt

# Define the file path where the plot will be saved
file_path = "/mnt/c/Users/User/Desktop/info/Compton/Simulation/plots/Interaction_probability/"

# Create a source object with specified photon energies (511 keV)
source = s.Source(energies={511: 1, 1274: 0})

number_of_photons = 100000
# Generate a list of test photons using the source
photons = source.testing_photons(number_of_photons)

# Function to simulate photons interacting with a target
def photon_in_target(photons, theta, number_of_steps=10):
    """
    Simulates photon propagation through a target and calculates the survival probability.
    
    :parameters photons: List of photon objects.
    :parameters theta: Incident angle of the photons in degrees.
    :parameters number_of_steps: Number of steps for photon propagation within the target.
    :returns: Probability of photon survival.
    """
    target_position = np.array([0, 3, 0])  
    target_width = 1 / np.abs(np.cos(np.deg2rad(theta))) 
    
    # Initialize target with material properties (Copper: Z=29, density=8.96 g/cm³, molar mass=63.55 g/mol)
    target = t.Target(position=target_position, radius=2, width=target_width, Z=29, density=8.96, molar_mass=63.55)

    distance = np.linalg.norm(target_position - source.position)
    step = target.width / number_of_steps

    n_dead_photons = 0 

    for photon in photons:
        # Propagate the photon towards the target
        photon = e.photon_propagation_to_target(photon, distance)

        # Simulate photon propagation inside the target in steps
        for j in range(number_of_steps):
            # Determine if the photon interacts with the target
            if random.uniform(0, 1) < i.interaction_probability(photon, step, target):
                n_dead_photons += 1  # Count photons that interact
                break  # Exit the loop if the photon scatters

    # Calculate the probability of photon survival
    Prob_survive = 1 - n_dead_photons / len(photons)

    return Prob_survive

thetas = np.arange(0, 80, 10)  
Probs_survive_compton = []  
Probs_survive = []  

for theta in thetas:
    # Compute probability of survival without scattering
    Probs_survive.append(photon_in_target(photons, theta))

    compton_photons = photons.copy()
    for photon in compton_photons:
        # Update photon energy after Compton scattering at given angle
        photon.energy = photon.compton_scattering(np.deg2rad(theta))

    # Compute probability of survival after Compton scattering
    Probs_survive_compton.append(photon_in_target(compton_photons, theta))

# Plot the probability behavior against angles
plt.figure(figsize=(8, 5))
plt.plot(thetas, Probs_survive, marker='o', linestyle='-', color='b', label="Photon Survival Probability")
plt.plot(thetas, Probs_survive_compton, marker='o', linestyle='-', color='r', label="Photon Survival Probability after Compton Scattering")
plt.xlabel(r'$\theta$ (degrees)')
plt.ylabel('Survival Probability')
plt.title('Photon Survival Probability vs Scattering Angle')
plt.legend()
plt.grid()
plt.savefig(file_path + "survivle_photons.png")
