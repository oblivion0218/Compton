import numpy as np
import random  
import target as t  
import experiments as e 
import particles as p  
import source as s 
import interactions as i
import matplotlib.pyplot as plt

file_path = "/mnt/c/Users/User/Desktop/info/Compton/Simulation/plots/Interaction_probability/"

# Create a source object with specified photon energies (511 keV and 1274 keV)
source = s.Source(energies={511: 1, 1274: 0})

# Define the number of photons to simulate
number_of_photons = 100000
# Generate a list of test photons using the source
photons = source.testing_photons(number_of_photons)

# Function to simulate photons interacting with a target
def photon_in_target(photons, theta, number_of_steps=10):

    target_position = np.array([0, 3, 0])  # Ensure it's an array for calculations
    target_width = 1 / np.abs(np.cos(np.deg2rad(theta)))
    
    target = t.Target(position=target_position, radius=2, width=target_width, Z=29, density=8.96, molar_mass=63.55)

    # Calculate the distance between the source and the target
    distance = np.linalg.norm(target_position - source.position)
    
    # Step size for photon propagation (target width divided by number of steps)
    step = target.width / number_of_steps

    n_dead_photons = 0

    # Loop over each photon in the list
    for photon in photons:
        # Propagate the photon towards the target
        photon = e.photon_propagation_to_target(photon, distance)

        # Simulate photon propagation inside the target in steps
        for j in range(number_of_steps):

            if random.uniform(0, 1) < i.interaction_probability(photon, step, target):
                n_dead_photons += 1
                break  # Exit the loop if the photon scatters

    # Calculate the probability of photon survival
    Prob_survive = 1 - n_dead_photons / len(photons)

    return Prob_survive

# Define angles and compute probabilities
thetas = np.arange(0, 80, 10)  # Angles in degrees

Probs_survive_compton = []
Probs_survive = []

for theta in thetas:
    Probs_survive.append(photon_in_target(photons, theta))

    compton_photons = photons.copy()
    for photon in compton_photons:
        photon.energy = photon.compton_scattering(np.deg2rad(theta))

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
