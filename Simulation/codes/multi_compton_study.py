import numpy as np
import random
import libraries.target as t
import libraries.experiments as e
import libraries.particles as p
import libraries.source as s
import libraries.interactions as i
import matplotlib.pyplot as plt

# Create a source object with specified photon energies (511 keV and 1274 keV)
source = s.Source(energies={511: 1, 1274: 0})

# Define the number of photons to simulate
number_of_photons = 10000

# Generate a list of test photons using the source
photons = source.testing_photons(number_of_photons)

# Set the position of the target in 3D space
target_position = [0, 5, 0]

# Create target objects for different materials (Pb, Al, and Cu) with specified properties
target_Pb = t.Target(position=target_position, radius=2, width=0.27, Z=82, density=11.34, molar_mass=207.2)
target_Al = t.Target(position=target_position, radius=2, width=3.17, Z=13, density=2.7, molar_mass=26.98)
target_Cu = t.Target(position=target_position, radius=2, width=0.93, Z=29, density=8.96, molar_mass=63.55)

# Function to simulate multiple Compton interactions of photons within a target
def multi_compton(photons, target, number_of_steps=10):
    # Calculate the distance between the source and the target
    distance = np.linalg.norm(target_position - source.position)
    
    # Calculate the step size for photon propagation inside the target
    step = target.width / number_of_steps

    # Counter to track photons that undergo two Compton scatterings
    n_photons_2_compton = 0
    
    # Loop through each photon in the list
    for photon in photons:
        # Propagate the photon from the source to the target
        photon = e.photon_propagation_to_target(photon, distance)

        # Simulate photon interactions within the target in steps
        for j in range(number_of_steps):
            # Check if the photon undergoes a Compton interaction
            if random.uniform(0, 1) < i.interaction_probability(photon, step, target) and i.which_interaction(photon, target.Z) == "compton":
                # Calculate scattering angles and update photon direction and energy
                theta = photon.compton_angle()
                phi = random.uniform(0, 1)
                x = np.sin(theta) * np.cos(phi)
                y = np.cos(theta)
                z = np.sin(theta) * np.sin(phi)

                # Update the photon's position and direction after the first Compton scattering
                photon.position = [0, j * step + distance, 0]
                photon.direction = [x, y, z]
                photon.energy = photon.compton_scattering(theta)

                # Check if the photon remains within the target after scattering
                while target.is_in_target(photon.position):
                    # If another Compton scattering occurs, increment the counter
                    if random.uniform(0, 1) < i.interaction_probability(photon, step, target) and i.which_interaction(photon, target.Z) == "compton":
                        n_photons_2_compton += 1
                        break
                    # Continue propagating the photon within the target
                    e.photon_propagation_to_target(photon, step)
                break

    return n_photons_2_compton

# Perform scattering simulations for each target material
photons_scattering_Pb = multi_compton(photons, target_Pb, number_of_steps=100)
photons_scattering_Al = multi_compton(photons, target_Al, number_of_steps=100)
photons_scattering_Cu = multi_compton(photons, target_Cu, number_of_steps=100)

# Calculate the probability of photon interaction in each target material
Prob_interaction_Pb = photons_scattering_Pb / len(photons)
Prob_interaction_Al = photons_scattering_Al / len(photons)
Prob_interaction_Cu = photons_scattering_Cu / len(photons)

# Print the results for each target material (Pb, Al, Cu)
print(f"Pb\tProbability that the photon interacts in the target: {Prob_interaction_Pb}\n")
print(f"Al\tProbability that the photon interacts in the target: {Prob_interaction_Al}\n")
print(f"Cu\tProbability that the photon interacts in the target: {Prob_interaction_Cu}\n")

# Define a range of target widths for further analysis
target_widths = np.linspace(0, 5, 50)

# Placeholder lists to store interaction probabilities for each material
prob_interaction_Pb = []
prob_interaction_Al = []
prob_interaction_Cu = []

# Loop over different target widths to calculate probabilities
for width in target_widths:
    # Update target widths for each material
    target_Pb.width = width
    target_Al.width = width
    target_Cu.width = width

    # Simulate photon interactions for the current target width
    photons_scattering_Pb = multi_compton(photons, target_Pb, number_of_steps=100)
    photons_scattering_Al = multi_compton(photons, target_Al, number_of_steps=100)
    photons_scattering_Cu = multi_compton(photons, target_Cu, number_of_steps=100)
    
    # Compute and store interaction probabilities
    prob_interaction_Pb.append(photons_scattering_Pb / len(photons))
    prob_interaction_Al.append(photons_scattering_Al / len(photons))
    prob_interaction_Cu.append(photons_scattering_Cu / len(photons))
    
# Plot Probability of Interaction vs Target Width
plt.figure(figsize=(12, 8))
plt.plot(target_widths, prob_interaction_Pb, label='Pb', color='blue')
plt.plot(target_widths, prob_interaction_Al, label='Al', color='green')
plt.plot(target_widths, prob_interaction_Cu, label='Cu', color='red')
plt.xlabel('Target Width (cm)')
plt.ylabel('Probability of 2 Compton Interaction')
plt.title('Probability of 2 Compton Interaction vs Target Width')
plt.legend()

# Save the plot as an image file
plt.tight_layout()
plt.savefig("multi_compton_probability.png")
