from lib import detector as d
from lib import particles as p
from lib import interactions
import numpy as np
import random
import os
import matplotlib.pyplot as plt
from tqdm import tqdm


file_path = "/mnt/c/Users/User/Desktop/info/Compton/Simulations/Python_simulation/interaction_point/"

r_e = 2.817e-13  # Classical electron radius in cm
m_e = 511  # Electron rest mass energy in keV
alpha = 1 / 137  # Fine-structure constant (dimensionless)
N_a = 6.022e23 # mol^(-1)

#For a Cu target
rho = 8.96 #g/cm^3
Z = 29 #number of electrons
MM = 63.546 #g/mol


def compton_scattering(angle: float) -> float:
    """
    Calculate the Compton scattering energy for a given angle.

    :param angle: Scattering angle in radians.
    :return: Compton scattering energy.
    """
    return 511 / (2 - np.cos(angle))  # Energy in keV


def attenuation_length(energy: float, Z: int) -> float:
    """
    Calculate the attenuation length for a photon traveling through a material.

    :param attenuation_factor: Attenuation factor (inverse mean free path).
    :return: Attenuation length in cm.
    """
    photon = p.Photon(energy, [0, 1, 0]) 
    target = d.Target(([0, 5, 0], [0, 6, 0]), 3) 
    cross_section = interactions.cross_section_photoelectric(photon, Z) + interactions.cross_section_compton(photon, Z)
    mu = interactions.attenuation_factor(cross_section, target)
    return 1 / mu  # cm


# Calculus of the probability factor 
def probability_factor(x: float, l: float, theta: float) -> float:
    """
    Calculate the probability factor for a photon interaction.

    :param x: Length before interaction point in cm.
    :param l: Length after interaction point in cm.
    :param theta: Angle of the photon in radians.
    :return: Probability factor.
    """
    E = 511  # Energy in keV
    lam = attenuation_length(E, Z)
    E_prime = compton_scattering(theta)  # Scattered energy in keV
    lam_prime = attenuation_length(E_prime, Z)

    return np.exp(-x / lam) * np.exp(-l / lam_prime)


target = d.Target(([0, 5, 0], [0, 6, 0]), 3)
interacion_points = []

with open(file_path + "data/IP_simulation_35.txt", 'r') as f:
    for i, line in enumerate(f):
        # Clean the line by removing brackets and splitting by commas or spaces
        clean_line = line.strip().replace('[', '').replace(']', '')
        
        # Try different splitting methods since we don't know the exact format
        if ',' in clean_line:
            # If commas are present, split by commas
            values = [float(val.strip()) for val in clean_line.split(',')]
        else:
            # Otherwise split by whitespace
            values = [float(val) for val in clean_line.split()]
        
            # Make sure we have 3 values (x, y, z)
        if len(values) == 3:
            interacion_points.append(values)
        else:
            print(f"Warning: Line {i+1} in file IP_simulation_35.txt does not contain exactly 3 values: {line}")

# INPUT points --> lenghts before IP 
x = []
source_position = [0, 0, 0]  # Source position

for ip in interacion_points:
    in_point = target.find_exit_point(np.array(ip), np.array(source_position))
    x.append(np.linalg.norm(in_point - np.array(ip)))

# OUTPUT points --> lenghts after IP
prob_factor = []
angles = [35, 40, 50, 60]

for angle in tqdm(angles, desc="Calculating probability factors"):
    l = []
    angle = np.radians(35) 
    detector_radius = 2.54  # cm
    detector_distance = 27.54  # cm
    theta_max = np.arctan(detector_radius / detector_distance) 

    skip_idx = []
    for i, ip in enumerate(interacion_points):
        theta = random.uniform(0, theta_max)
        phi = random.uniform(0, 2 * np.pi)
        r = detector_distance * np.tan(theta)

        external_point = np.array([r * np.sin(phi), detector_distance, r * np.cos(phi)])
        rotation_center = np.array([0, 5.5, 0])  # Center of the target
        rotation_matrix = np.array([
                    [np.cos(angle), -np.sin(angle), 0],
                    [np.sin(angle), np.cos(angle), 0],
                    [0, 0, 1]
                ])
        rotated_external_point = np.dot(rotation_matrix, external_point - rotation_center) + rotation_center

        out_point = target.find_exit_point(np.array(ip), np.array(rotated_external_point))

        if out_point is None:
            skip_idx.append(i)
            continue

        else: 
            l.append(np.linalg.norm(out_point - np.array(ip)))


    interacion_points_def = [ip for i, ip in enumerate(interacion_points) if i not in skip_idx]
    x_def = [x_val for i, x_val in enumerate(x) if i not in skip_idx]
    l_def = l

    correction = 0
    for l, x in zip(l_def, x_def):
        correction += probability_factor(x, l, angle)

    correction /= len(l_def)
    prob_factor.append(correction)

print(prob_factor)