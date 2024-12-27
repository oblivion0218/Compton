import numpy as np
import matplotlib.pyplot as plt
import random
from scipy.integrate import quad

# Constants and initial parameters
#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

source_position = np.array([0, 0, 8])  # Position of the source

detector_params = {
    "ugo": {"radius": 1.27, "width": 2.54, "position": np.array([0, 5, 8])},
    "franco": {"radius": 2.54, "width": 5.08, "position": np.array([0, -5, 8])},
}

density = 3.67  # g/cm^3 (Density of NaI, detector material)
Z = 47  # Approximate atomic number of NaI (mean of Na and I)

# Physical constants
r_e = 2.817e-15  # Classical electron radius in m
m_e = 511  # Rest mass energy of the electron in keV
alpha = 1 / 137  # Fine-structure constant

# Class Photon
#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

class Photon:
    def __init__(self, energy: float, direction: list[float]):
        self.energy = energy
        self.direction = np.array(direction)

    def __str__(self):
        return f"Photon with energy {self.energy} keV and direction {self.direction}"

    def __repr__(self):
        return f"Photon({self.energy}, {self.direction})"

    def propagation(self, distance: float, source_pos: np.ndarray = source_position) -> np.ndarray:
        """Calculate the position of the photon after traveling a given distance."""
        return source_pos + self.direction * distance

    def is_in_detector(self, detector_params: dict, distance: float) -> bool:
        """Check if the photon hits a specific detector."""
        detector_pos = detector_params["position"]
        detector_radius = detector_params["radius"]

        hit_point = self.propagation(distance)
        within_radius = np.linalg.norm(hit_point[[0, 2]] - detector_pos[[0, 2]]) <= detector_radius
        same_side = np.sign(hit_point[1]) == np.sign(detector_pos[1])

        return within_radius and same_side

    def interaction_probability(self, detector_params: dict, distance: float) -> float:
        """Calculate the interaction probability within the detector."""
        if self.is_in_detector(detector_params, distance):
            mu = 0.2 if self.energy <= 511 else 0.1  # Linear attenuation coefficient in cm^-1
            return 1 - np.exp(-distance * mu)
        return 0

# Functions to initialize photons    
#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

def photons(number_of_photons: int) -> list:
    """Generate a list of photons with random directions and energies."""
    phi = np.random.uniform(0, 2 * np.pi, number_of_photons)
    theta = np.random.uniform(0, np.pi, number_of_photons)

    directions = np.vstack([
        np.sin(theta) * np.cos(phi),
        np.sin(theta) * np.sin(phi),
        np.cos(theta)
    ]).T

    energies = np.random.choice([511, 1274], size=number_of_photons, p=[0.903, 0.097])
    return [Photon(energy, direction) for energy, direction in zip(energies, directions)]

# Cross-sections
#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*



#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
# Main Simulation   ---> Ovviamente non ha senso 
#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

# def simulate_photons(number_of_photons: int, distance: float):
#     """Simulate photons and check interactions with detectors."""
#     photon_list = photons(number_of_photons)

#     for photon in photon_list:
#         if photon.is_in_detector(detector_params["ugo"], distance):
#             interaction_prob = photon.interaction_probability(detector_params["ugo"], distance)
#             print(f"{photon} hits the Ugo detector with interaction probability {interaction_prob:.2f}.")

#         if photon.is_in_detector(detector_params["franco"], distance):
#             interaction_prob = photon.interaction_probability(detector_params["franco"], distance)
#             print(f"{photon} hits the Franco detector with interaction probability {interaction_prob:.2f}.")

# # Run the simulation
# simulate_photons(1000, 10)
