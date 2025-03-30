import numpy as np
import random  
from lib import detector as d
from lib import experiments as e 
from lib import particles as p  
from lib import source as s 
from lib import interactions as i
import matplotlib.pyplot as plt
from tqdm import tqdm

file_path = "/mnt/c/Users/User/Desktop/info/Compton/Simulation/plots/Interaction_probability/"

# Create a source object with specified photon energies (511 keV and 1274 keV)
source = s.Source(energies={511: 1, 1274: 0})

# Define the number of photons to simulate
number_of_photons = 100000
# Generate a list of test photons using the source
photons = source.testing_photons(number_of_photons)

# Set the position of the target in 3D space
target_position = ([0, 5, 0], [0, 10.08, 0])  # Example position (arbitrary values)

# Create target objects for different materials (Pb, Al, and Cu) with specified properties
target_Pb = d.Target(position=target_position, radius=2.54, Z=82, density=11.34, molar_mass=207.2)
target_Al = d.Target(position=target_position, radius=2.54, Z=13, density=2.7, molar_mass=26.98)
target_Cu = d.Target(position=target_position, radius=2.54, Z=29, density=8.96, molar_mass=63.55)


# Function to simulate photons interacting with a target
def photon_in_target(photons, target, number_of_steps=10):
    # Calculate the distance between the source and the target
    distance = np.linalg.norm(target_position - source.position)
    
    # Step size for photon propagation (target width divided by number of steps)
    step = target.width / number_of_steps

    n_scattering_photons = 0

    # Loop over each photon in the list
    for photon in photons:
        # Propagate the photon towards the target
        photon = e.photon_propagation_to_target(photon, distance)

        # Simulate photon propagation inside the target in steps
        for j in range(number_of_steps):
            # If a random number is less than the interaction probability, the photon scatters
            if random.uniform(0, 1) < i.interaction_probability(photon, step, target):
                n_scattering_photons += 1
                break  # Exit the loop if the photon scatters

    return n_scattering_photons

target_widths = np.linspace(0, 2, 50)

# Placeholder lists to store probabilities for each material
prob_interaction_Pb = []
prob_interaction_Al = []
prob_interaction_Cu = []

prob_compton_Pb = []
prob_compton_Al = []
prob_compton_Cu = []

# Loop over different target widths and calculate probabilities
for width in tqdm(target_widths, desc="Calculating probabilities", unit="width"):
    target_Pb.width = width
    target_Al.width = width
    target_Cu.width = width

    # Simulate interactions for current width
    photons_scattering_Pb = photon_in_target(photons, target_Pb, number_of_steps=100)
    photons_scattering_Al = photon_in_target(photons, target_Al, number_of_steps=100)
    photons_scattering_Cu = photon_in_target(photons, target_Cu, number_of_steps=100)
    
    # Compute interaction probabilities
    prob_interaction_Pb.append(photons_scattering_Pb / len(photons))
    prob_interaction_Al.append(photons_scattering_Al / len(photons))
    prob_interaction_Cu.append(photons_scattering_Cu / len(photons))
    
    # Compute Compton probabilities (assuming cross_section functions are correctly implemented)
    prob_compton_Pb.append(
        prob_interaction_Pb[-1] * i.cross_section_compton(photons[0], target_Pb.Z) /
        (i.cross_section_compton(photons[0], target_Pb.Z) + i.cross_section_photoelectric(photons[0], target_Pb.Z))
    )
    prob_compton_Al.append(
        prob_interaction_Al[-1] * i.cross_section_compton(photons[0], target_Al.Z) /
        (i.cross_section_compton(photons[0], target_Al.Z) + i.cross_section_photoelectric(photons[0], target_Al.Z))
    )
    prob_compton_Cu.append(
        prob_interaction_Cu[-1] * i.cross_section_compton(photons[0], target_Cu.Z) /
        (i.cross_section_compton(photons[0], target_Cu.Z) + i.cross_section_photoelectric(photons[0], target_Cu.Z))
    )

# Plot Probability of Interaction vs Target Width
plt.figure(figsize=(12, 5))
plt.subplot(1, 2, 1)
plt.plot(target_widths, prob_interaction_Pb, label='Pb', color='blue')
plt.plot(target_widths, prob_interaction_Al, label='Al', color='green')
plt.plot(target_widths, prob_interaction_Cu, label='Cu', color='red')
plt.xlabel('Target Width (cm)')
plt.ylabel('Probability of Interaction')
plt.title('Probability of Interaction vs Target Width')
plt.grid(True)
plt.legend()

# Plot Probability of Compton Interaction vs Target Width
plt.subplot(1, 2, 2)
plt.plot(target_widths, prob_compton_Pb, label='Pb', color='blue')
plt.plot(target_widths, prob_compton_Al, label='Al', color='green')
plt.plot(target_widths, prob_compton_Cu, label='Cu', color='red')
plt.xlabel('Target Width (cm)')
plt.ylabel('Probability of Compton Interaction')
plt.title('Probability of Compton Interaction vs Target Width')
plt.grid(True)
plt.legend()

plt.tight_layout()
plt.savefig("interaction_probability_02.png")