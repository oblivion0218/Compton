import numpy as np
import matplotlib.pyplot as plt
import random
from scipy.integrate import quad


# Constants and initial parameters
#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

source_position = np.array([0, 0, 8])  # Position of the source (0,0,0)

# Detector parameters
detector_ugo_radius = 1.27  # Radius in cm (1 inch = 2.54 cm)
detector_ugo_width = 2.54  # Width in cm
detector_franco_radius = 2.54  # Radius in cm
detector_franco_width = 5.08  # Width in cm

ugo_position = np.array([0, 5, 8])  # Distance from source (5 cm along y-axis)
franco_position = np.array([0, -5, 8])  # Distance from source (-5 cm along y-axis)

density = 3.67  # g/cm^3 (Density of NaI, detector material)
Z = 47  # Atomic number of NaI (detector material), it is a sort of mean between Z = 11 (Na) and Z = 53 (I)

# Physical constants
r_e = 2.817e-15  # Classical electron radius in m
m_e = 511  # Rest mass energy of the electron in keV
alpha = 1 / 137  # Fine-structure constant

#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

class Photon:
    def __init__(self, energy, direction):
        self.energy = energy
        self.direction = direction
    
    def __str__(self):
        return f"Photon with energy {self.energy} keV and direction {self.direction}"

    def __repr__(self):
        return f"Photon({self.energy}, {self.direction})"

    def propagation(self, distance, source_pos=source_position): # Aggiungi la riflessione dei fotoni che hanno z < 0
        x = source_pos[0] + self.direction[0] * distance
        y = source_pos[1] + self.direction[1] * distance
        z = source_pos[2] + self.direction[2] * distance
        return np.array([x, y, z])

    def is_in_detector(self, detector, distance, source_pos=source_position):
        if detector == "ugo":
            detector_pos = ugo_position
            detector_radius = detector_ugo_radius
        elif detector == "franco":
            detector_pos = franco_position
            detector_radius = detector_franco_radius

        hit_point = self.propagation(distance)
        if detector_pos[1] * hit_point[1] > 0 and np.sqrt((hit_point[0] - source_pos[0]) ** 2 + (hit_point[2] - source_pos[2]) ** 2) <= detector_radius:
            return True
        else:
            return False
        
    def interaction_probability(self, detector, distance):
        if self.is_in_detector(detector, distance):
            mu = 0

            if self.energy == 511:
                mu = 0.2 # cm^-1    
            elif self.energy == 1274:
                mu = 0.1 # cm^-1
            
            return (1 - np.exp(- distance * mu))
        
        else:
            return 0
        
# Functions to inizialized the photons    
#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

# Function that chooses the type of photon to be emitted
def which_photon():
    photon_type = random.uniform(0, 100)
    if photon_type <= 90.3:
        return 511
    else:
        return 1274

# Function to generate a list of photons object with their energy and direction
def photons(number_of_photons):
    photon_list = []  
    for _ in range(number_of_photons):
        phi = np.random.uniform(0, 2 * np.pi)
        theta = np.random.uniform(0, np.pi)

        x = np.sin(theta) * np.cos(phi)
        y = np.sin(theta) * np.sin(phi)
        z = np.cos(theta)

        energy = which_photon()
        direction = [x, y, z]
        photon_list.append(Photon(energy, direction))  

    return photon_list 

# Cross-sections
#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

