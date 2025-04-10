import numpy as np
import matplotlib.pyplot as plt
import random
from tqdm import tqdm
from lib import interactions as i
from lib import detector as d
from lib import particles as p
from lib import experiments as e
from lib import source as s
from lib import visualization as v

# File path to save the output spectrum plot
file_path = "/mnt/c/Users/User/Desktop/info/Gamma-simulation/simulated_events/"      #ANDRE

# Object initialization
# *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
source = s.Source({511: 1, 1274: 0})  # Create a source object

# 110°
target = d.Target(([0, 5, 0], [0, 6, 0]), 3)  # Create a target object
target.rotate((-7/36) * np.pi, [0, 5.5, 0], "z")  # Rotate the target

detector = d.Detector(([0, 30.5, 0], [0, 35.58, 0]), 2.54, 0.0695)  # Detector "Franco"
detector.rotate((11/18) * np.pi, [0, 5.5, 0], "z")  # Rotate the detector 

# Initial parameters
# *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
step = 0.1
N_cycles = 20
number_of_photons = 1000000

r_gate = 1.27
d_gate_source = 16
theta_gate = np.arctan(r_gate/d_gate_source)

def new_direction(photon):
    """
    Calculate the new direction of a photon after a Compton scattering interaction.

    :param photon: Photon object.
    :return: New direction vector and scattering angle.
    """
    theta = photon.compton_angle()
    phi = np.random.uniform(0, 2 * np.pi)
    x = np.sin(theta) * np.cos(phi)
    y = np.cos(theta)
    z = np.sin(theta) * np.sin(phi)
    return np.array([x, y, z]), theta

for j in tqdm(range(N_cycles), desc="Simulating cycles", unit="cycle"):
    # Simulation of photons emitted from the source 
    # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    photons = source.photon_emission(number_of_photons, theta_gate, 2 * np.pi, axis="y", forward_backward=False) 

    [e.photon_propagation_to_target(photon, target) for photon in photons] 
    # v.visualization_3D_plotly("photons.html", [detector], photons, source, target)
    [photon.propagation(step) for photon in photons]  # Propagate the photons

    photons_out_of_target = []
    n_compton_interactions = [0, 0, 0, 0]  # [first compton, second compton, third compton, fourth compton]
    n_photoelectric_interactions = [0, 0, 0, 0]  # [first photoelectric, second photoelectric, third photoelectric, fourth photoelectric]
    n_interaction = [0, 0, 0, 0]  # [first interaction, second interaction, third interaction, fourth interaction]

    # Photons interacting with the target
    # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    for photon in tqdm(photons, desc="Simulating photons", unit="photon"):
        while target.is_inside(photon.position):
            # I interaction
            if i.interaction_probability(photon, step, target) > random.uniform(0, 1):
                n_interaction[0] += 1
                # I photoelectric interaction
                if i.which_interaction(photon, target.Z) == "photoelectric":
                    n_photoelectric_interactions[0] += 1
                    photon.energy = 0 
                    break
                # I Compton interaction
                else:
                    n_compton_interactions[0] += 1
                    new_dir, scattering_angle = new_direction(photon)
                    photon.direction = new_dir
                    photon.energy = photon.compton_scattering(scattering_angle)

                    photon.propagation(step)

                    while target.is_inside(photon.position):
                        # II interaction
                        if i.interaction_probability(photon, step, target) > random.uniform(0, 1):
                            n_interaction[1] += 1
                            # II photoelectric interaction
                            if i.which_interaction(photon, target.Z) == "photoelectric":
                                n_photoelectric_interactions[1] += 1
                                photon.energy = 0 
                                break
                            # II Compton interaction
                            else:
                                n_compton_interactions[1] += 1
                                new_dir, scattering_angle = new_direction(photon)
                                photon.direction = new_dir
                                photon.energy = photon.compton_scattering(scattering_angle)

                                photon.propagation(step)

                                while target.is_inside(photon.position):
                                    # III interaction
                                    if i.interaction_probability(photon, step, target) > random.uniform(0, 1):
                                        n_interaction[2] += 1
                                        # III photoelectric interaction
                                        if i.which_interaction(photon, target.Z) == "photoelectric":
                                            n_photoelectric_interactions[2] += 1
                                            photon.energy = 0 
                                            break
                                        # III Compton interaction
                                        else:
                                            n_compton_interactions[2] += 1
                                            new_dir, scattering_angle = new_direction(photon)
                                            photon.direction = new_dir
                                            photon.energy = photon.compton_scattering(scattering_angle)

                                            photon.propagation(step)

                                            while target.is_inside(photon.position):
                                                # IV interaction
                                                if i.interaction_probability(photon, step, target) > random.uniform(0, 1):
                                                    n_interaction[3] += 1
                                                    # IV photoelectric interaction
                                                    if i.which_interaction(photon, target.Z) == "photoelectric":
                                                        n_photoelectric_interactions[3] += 1
                                                        photon.energy = 0 
                                                        break
                                                    # IV Compton interaction
                                                    else:
                                                        n_compton_interactions[3] += 1
                                                        new_dir, scattering_angle = new_direction(photon)
                                                        photon.direction = new_dir
                                                        photon.energy = photon.compton_scattering(scattering_angle)
                                                        break
                                                # close IV interaction
                                                else:
                                                    photon.propagation(step)
                                            break
                                    # close III interaction
                                    else:
                                        photon.propagation(step)  
                                break
                        # close II interaction
                        else:
                            photon.propagation(step) 
                    break            
            # close I interaction
            else:
                photon.propagation(step)

        if photon.energy:
            photons_out_of_target.append(photon)    

    # Print the results on screen
    # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    print(f"\nNumber of photons that left the target: {len(photons_out_of_target)}\n")

    print(f"Number of photons that interacted wi8th the target: {n_interaction[0]}")
    print(f"Number of photons that interacted with the target (2nd): {n_interaction[1]}")
    print(f"Number of photons that interacted with the target (3rd): {n_interaction[2]}")
    print(f"Number of photons that interacted with the target (4th): {n_interaction[3]}\n")

    print(f"Number of photoelectric interactions: {n_photoelectric_interactions[0]}")
    print(f"Number of photoelectric interactions (2nd): {n_photoelectric_interactions[1]}")
    print(f"Number of photoelectric interactions (3rd): {n_photoelectric_interactions[2]}")
    print(f"Number of photoelectric interactions (4th): {n_photoelectric_interactions[3]}\n")

    print(f"Number of Compton interactions: {n_compton_interactions[0]}")
    print(f"Number of Compton interactions (2nd): {n_compton_interactions[1]}")
    print(f"Number of Compton interactions (3rd): {n_compton_interactions[2]}")
    print(f"Number of Compton interactions (4th): {n_compton_interactions[3]}\n")

    # Photons that reached the detector
    # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    detector_axis = detector.principal_axis()/np.linalg.norm(detector.principal_axis())
    theta_detector = np.arctan(detector.radius / np.linalg.norm([0, 5.5, 0] - detector.position[0]))

    photons_to_detector = []

    for photon in photons_out_of_target:
        distance = np.linalg.norm(photon.position - detector.position[0]) # distance to the detector
        # Check if the photon is within the detector's acceptance angle
        if np.arccos(np.dot(photon.direction, detector_axis)) < theta_detector:
            e.photon_propagation_to_target(photon, detector)
            photons_to_detector.append(photon)
        else:
            photon.propagation(distance)

    # v.visualization_3D_plotly(file_path + "3D_visualization/survival_photons.html", [detector], photons_out_of_target, source, target)
    print(f"Number of photons that reached the detector: {len(photons_to_detector)}")
    # v.visualization_3D_plotly(file_path + "3D_visualization/photons_to_detector.html", [detector], photons_to_detector, source, target)

    # Photons observed by the detector
    # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    detected_energies = [] #array for eletron's energy as detected by the detector
    true_detected_energies = [] #array for electron's real energy 
    for photon in photons_to_detector:
        energy = e.gamma_detection(photon, detector, step=step, true_energy=True)
        if energy > 0:
            detected_energies.append(detector.resolution(energy))
            true_detected_energies.append(energy)

    print(f"Number of photons observed by the detector: {len(detected_energies)}\n")

    # Print the results in a txt file
    # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    with open(f"{j}_detected_energies.txt", "w") as f:
        f.write(f"Number of photons that left the target: {len(photons_out_of_target)}\n")

        f.write(f"Number of photons that interacted with the target: {n_interaction[0]}\n") 
        f.write(f"Number of photons that interacted with the target (2nd): {n_interaction[1]}\n")
        f.write(f"Number of photons that interacted with the target (3rd): {n_interaction[2]}\n")
        f.write(f"Number of photons that interacted with the target (4th): {n_interaction[3]}\n")

        f.write(f"Number of photoelectric interactions: {n_photoelectric_interactions[0]}\n")
        f.write(f"Number of photoelectric interactions (2nd): {n_photoelectric_interactions[1]}\n")
        f.write(f"Number of photoelectric interactions (3rd): {n_photoelectric_interactions[2]}\n")
        f.write(f"Number of photoelectric interactions (4th): {n_photoelectric_interactions[3]}\n")

        f.write(f"Number of Compton interactions: {n_compton_interactions[0]}\n")
        f.write(f"Number of Compton interactions (2nd): {n_compton_interactions[1]}\n")
        f.write(f"Number of Compton interactions (3rd): {n_compton_interactions[2]}\n")
        f.write(f"Number of Compton interactions (4th): {n_compton_interactions[3]}\n")

        f.write(f"Number of photons that reached the detector: {len(photons_to_detector)}\n")
        f.write(f"Number of photons observed by the detector: {len(detected_energies)}\n")

        f.write("Detected Energies:\n")
        for energy in detected_energies:
            f.write(f"{energy}\n")
        f.write("True Detected Energies:\n")
        for energy in true_detected_energies:
            f.write(f"{energy}\n")


# Simulation parameters for 511 keV
# *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
p1 = p.Photon(511, [0, 0, 0])
sigma_photoelectric = i.cross_section_photoelectric(p1, target.Z)
sigma_compton = i.cross_section_compton(p1, target.Z)
sigma = sigma_photoelectric + sigma_compton
mu = i.attenuation_factor(sigma, target)
prob = i.interaction_probability(p1, 0.1, target)
lamda = 1 / mu
print("SIMULATION PARAMETERS")
print(f"Total cross section: {sigma}")
print(f"Photoelectric cross section: {sigma_photoelectric}")
print(f"Compton cross section: {sigma_compton}")
print(f"Attenuation factor: {mu}")
print(f"Mean free path: {lamda}")
print(f"Interaction probability: {prob}")


# Results for 1,000,000 photons
# *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
"""
Number of photons that left the target: 825220

Number of photons that interacted wi8th the target: 635876
Number of photons that interacted with the target (2nd): 326167
Number of photons that interacted with the target (3rd): 157921
Number of photons that interacted with the target (4th): 67387

Number of photoelectric interactions: 29424
Number of photoelectric interactions (2nd): 57407
Number of photoelectric interactions (3rd): 54453
Number of photoelectric interactions (4th): 33496

Number of Compton interactions: 606452
Number of Compton interactions (2nd): 268760
Number of Compton interactions (3rd): 103468
Number of Compton interactions (4th): 33891

Number of photons that reached the detector: 442
Number of photons observed by the detector: 328


SIMULATION PARAMETERS for 511 keV
Total cross section: 8.706368970239455e-24 cm^2
Photoelectric cross section: 4.022561295784187e-25 cm^2
Compton cross section: 8.304112840661037e-24 cm^2
Attenuation factor: 0.7371980163078985 cm^-1
Mean free path: 1.3564876435890185 cm
Interaction probability: 0.07106805740563149
"""