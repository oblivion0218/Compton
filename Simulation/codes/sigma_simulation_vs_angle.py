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

# Initialize detectors with their respective positions and dimensions
detector = d.Detector(([0, 30.5, 0], [0, 35.58, 0]), 2.54, 0.0695)  # Detector "Franco"

number_of_photons = 1000000  # Number of photons to simulate

# Simulate photon emission from a source
source = s.Source({511: 1, 1274: 0})  # Create a source object

detector.rotate((1/2) * np.pi, [0, 5.5, 0], "z")

target = d.Target(([0, 5, 0], [0, 6, 0]), 3)  # Create a target object
target.rotate((-1/4) * np.pi, [0, 5.5, 0], "z")  # Rotate the target

step = 0.1

for j in tqdm(range(1), desc="N_cycles", unit="iteration"): 
    photons = source.photon_emission(number_of_photons, np.arctan(1.27/16), 2 * np.pi, axis="y", forward_backward=False)  # Simulate photon emission
    [e.photon_propagation_to_target(photon, target) for photon in photons]  # Propagate the photons
    v.visualization_3D_plotly(file_path + "3D_visualization/photons.html", [detector], photons, source, target)

    photons_out_of_target = []
    n_interacting_photons = 0
    n_compton_interactions = [0, 0, 0]  # [first compton, second compton]
    n_photoelectric_interactions = 0

    def new_direction(photon):
        theta = photon.compton_angle()
        phi = np.random.uniform(0, 2 * np.pi)
        x = np.sin(theta) * np.cos(phi)
        y = np.cos(theta)
        z = np.sin(theta) * np.sin(phi)
        return np.array([x, y, z]), theta

    # Photon propagation inside the target
    for photon in tqdm(photons, desc="Simulating photons", unit="photon"):
        photon.propagation(step)

        while target.is_inside(photon.position):
            if i.interaction_probability(photon, step, target) > random.uniform(0, 1):
                n_interacting_photons += 1
                if i.which_interaction(photon, target.Z) == "compton":
                    n_compton_interactions[0] += 1

                    photon.direction, theta1 = new_direction(photon)
                    photon.energy = photon.compton_scattering(theta1)

                    photon.propagation(step)

                    # Second compton scattering
                    while target.is_inside(photon.position):
                        if i.interaction_probability(photon, step, target) > random.uniform(0, 1):
                            if i.which_interaction(photon, target.Z) == "compton":
                                n_compton_interactions[1] += 1

                                photon.direction, theta2 = new_direction(photon)
                                photon.energy = photon.compton_scattering(theta2)

                                photon.propagation(step)

                                # Third compton scattering
                                while target.is_inside(photon.position):   
                                    if i.interaction_probability(photon, step, target) > random.uniform(0, 1):
                                        if i.which_interaction(photon, target.Z) == "compton":
                                            n_compton_interactions[2] += 1

                                            photon.direction, theta3 = new_direction(photon)
                                            photon.energy = photon.compton_scattering(theta3)
                                            break
                                        else:
                                            photon.energy = 0
                                            break
                                    break
                            else:
                                photon.energy = 0
                                break 
                            break             
                        else:
                            break
                    break
                else: 
                    n_photoelectric_interactions += 1
                    photon.energy = 0
                    break    
        if photon.energy > 0:
            photons_out_of_target.append(photon)

    print(f"\nNumber of photons that left the target: {len(photons_out_of_target)}")
    print(f"Number of photons that interacted with the target: {n_interacting_photons}")
    print(f"Number of photoelectric interactions: {n_photoelectric_interactions}")
    print(f"Number of Compton interactions: {n_compton_interactions[0]}")
    print(f"Number of Compton interactions (2nd): {n_compton_interactions[1]}")
    print(f"Number of Compton interactions (3rd): {n_compton_interactions[2]}\n")

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

    v.visualization_3D_plotly(file_path + "3D_visualization/survival_photons.html", [detector], photons_out_of_target, source, target)
    print(f"Number of photons that reached the detector: {len(photons_to_detector)}")
    v.visualization_3D_plotly(file_path + "3D_visualization/photons_to_detector.html", [detector], photons_to_detector, source, target)

    detected_energies = [] #array for eletron's energy as detected by the detector
    true_detected_energies = [] #array for electron's real energy 
    for photon in photons_to_detector:
        energy = e.gamma_detection(photon, detector, step=step, true_energy=True)
        if energy > 0:
            detected_energies.append(detector.resolution(energy))
            true_detected_energies.append(energy)

    # v.plot_energy_spectrum(detected_energies, "Energies.png")
    print(f"Number of photons observed by the detector: {len(detected_energies)}")

    with open(f"{j}_detected_energies.txt", "w") as f:
        f.write(f"Number of photons that left the target: {len(photons_out_of_target)}\n")
        f.write(f"Number of photons that interacted with the target: {n_interacting_photons}\n")
        f.write(f"Number of photoelectric interactions: {n_photoelectric_interactions}\n")
        f.write(f"Number of Compton interactions: {n_compton_interactions[0]}\n")
        f.write(f"Number of Compton interactions (2nd): {n_compton_interactions[1]}\n")
        f.write(f"Number of Compton interactions (3rd): {n_compton_interactions[2]}\n")
        f.write(f"Number of photons that reached the detector: {len(photons_to_detector)}\n")
        f.write(f"Number of photons observed by the detector: {len(detected_energies)}\n")
        f.write("Detected Energies:\n")
        for energy in detected_energies:
            f.write(f"{energy}\n")
        f.write("True Detected Energies:\n")
        for energy in true_detected_energies:
            f.write(f"{energy}\n")

# Results for 1,000,000 photons
"""
Number of photons that left the target: 963976
Number of photons that interacted with the target: 523522
Number of photoelectric interactions: 24046
Number of Compton interactions: 499476
Number of Compton interactions (2nd): 36435
Number of Compton interactions (3rd): 2916

Number of photons that reached the detector: 1179
Number of photons observed by the detector: 1179
"""

