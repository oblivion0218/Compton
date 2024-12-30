import numpy as np
import matplotlib.pyplot as plt
import random
from scipy.integrate import quad

# Constants and initial parameters
#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

source_position = np.array([0, 0, 0])  # Position of the source

detector_params = {
    "ugo": {"radius": 1.27, "width": 2.54, "position": np.array([0, 5, 0]), "resolution": 0.0695},
    "franco": {"radius": 2.54, "width": 5.08, "position": np.array([0, -5, 0]), "resolution": 0.0903},
}

density = 3.67  # g/cm^3 (Density of NaI, detector material)
Z = 47  # Approximate atomic number of NaI (mean of Na and I)

# Physical constants
bn = 1e-28  # Barn in m^2
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
        return self.direction * distance

    def will_be_in_detector(self, detector_params: dict, source_position: np.ndarray = source_position) -> bool:
        """Check if the photon will hit a specific detector."""

        detector_pos = detector_params["position"]
        detector_radius = detector_params["radius"]

        distance = np.linalg.norm(detector_pos - source_position)
        detector_theta = np.arctan(detector_radius / distance)

        hit_point = self.propagation(distance)
        radius_hit = np.linalg.norm(hit_point[[0, 2]] - detector_pos[[0, 2]])
        theta_hit = np.arctan(radius_hit / hit_point[1])

        within_radius = theta_hit <= detector_theta
        same_side = np.sign(hit_point[1]) == np.sign(detector_pos[1])

        return within_radius and same_side
    
    def is_in_detector(self, detector_params: dict, distance: float) -> bool: # Controlla
        """Check if the photon hits a specific detector."""
        detector_pos = detector_params["position"]
        detector_radius = detector_params["radius"]

        hit_point = self.propagation(distance)
        within_radius = np.linalg.norm(hit_point[[0, 2]] - detector_pos[[0, 2]]) <= detector_radius
        within_width = hit_point[1] >= detector_pos[1] and hit_point[1] <= detector_pos[1] + detector_params["width"]

        return within_radius and within_width

    def compton_scattering(self, angle: float) -> float:
        """Calculate the energy of the photon after Compton scattering."""
        return self.energy / (1 + (self.energy / m_e) * (1 - np.cos(angle)))        

    def klein_nishina(self, angle: float) -> float:
        """Calculate the Klein-Nishina differential cross-section."""
        r = self.compton_scattering(angle) / self.energy
        c = alpha ** 2 / (2 * m_e ** 2)
        return c * r ** 2 * (r + 1 / r - np.sin(angle) ** 2)
    
def pdf_compton(photon: Photon, angle: float) -> float:
    """Calculate the probability density function for Compton scattering."""
    return photon.klein_nishina(angle) / quad(photon.klein_nishina, 0, np.pi)[0]

def random_angles(photon: Photon) -> float:
    """Generate random angles using rejection sampling."""
    max_pdf = 0.9  # Adjust as needed; this is M for scaling
    while True:
        theta = np.random.uniform(0, np.pi)
        u = np.random.uniform(0, max_pdf)
        # Accept or reject
        if u <= pdf_compton(photon, theta):
            return theta


# Functions to initialize photons    
#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

def photons(number_of_photons) -> list:
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

def testing_photons(number_of_photons) -> list:
    """Generate a list of photons with random directions and energies."""
    
    directions_po = []
    for i in range(number_of_photons):
        directions_po.append(np.array([0, 1, 0]))
    directions = np.array(directions_po)
    energies = np.random.choice([511, 1274], size=number_of_photons, p=[0.903, 0.097])
    return [Photon(energy, direction) for energy, direction in zip(energies, directions)]

# Cross-sections form "A Modern Primer to Particle and Nuclear Physycs" by F.Terranova page 69 - 71
# I don't consider the pair production because it's not relevant for the energy range we are considering
#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

def cross_section_photoelectric(energy: float, detector_Z: int = Z) -> float: # It's well defined for energy << m_e
    """Calculate the photoelectric cross-section for a given energy in the detector."""
    epsilon = energy / m_e
    if epsilon <= 1: # It's well defined for energy << m_e
        c = 0.665 * np.sqrt(32) * alpha ** 4 * bn
        # c = (16 / 3) * np.sqrt(2) * np.pi * alpha ** 4 * r_e ** 2 # wikipedia
        return c * detector_Z ** 5 / (epsilon) ** 3.5
    else: # It's well defined for energy >> m_e
        c = 0.665 * (3 / 2) * alpha ** 4 * bn
        return c * detector_Z ** 5 / (epsilon)
    # we are studing a borderline case so the simulation is not very accurate    

def cross_section_compton(energy: float, detector_Z: int = Z) -> float:
    """Calculate the Compton cross-section for a given energy."""
    cross_section_thompson = (8 / 3) * np.pi * r_e ** 2
    epsilon = energy / m_e
    c = (3 / 4) * cross_section_thompson
    return c * detector_Z * (((1 + epsilon) / epsilon ** 2) * ((2 *(1 + epsilon) / (1 + 2 * epsilon)) - (np.log(1 + 2 * epsilon) / epsilon)) + (np.log(1 + 2 * epsilon) / (2 * epsilon)) - (1 + 3 * epsilon) / (1 + 2 * epsilon) ** 2)

#-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
# 511: σ_pe = 2.4490838560315705e-28 ; σ_com = 1.3458389776243748e-27
# σ_pe / σ_tot = 0.15395806478424864 ; σ_com / σ_tot = 0.8460419352157514 (σ_tot = σ_pe + σ_com)
#-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
# 1274: σ_pe = 1.0008841411785894e-29 ; σ_com = 8.782021780069243e-28
# σ_pe / σ_tot = 0.01126854001241302 ; σ_com / σ_tot = 0.9887314599875869(σ_tot = σ_pe + σ_com)
#-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

def interaction_probability(photon: Photon, distance: float) -> float:
    """Calculate the interaction probability within the detector."""

    attenuation_data = { # cm^-1
    10: 450,
    50: 25.0,
    100: 11.0,
    200: 4.5,
    511: 2.1,
    1000: 1.0,
    1274: 0.8
    }

    for attenuation_energy in attenuation_data:
        if photon.energy <= attenuation_energy:
            mu = attenuation_data[attenuation_energy]
            return 1 - np.exp(- mu * distance)

print(interaction_probability(Photon(511, [0, 1, 0]), 0.1))
# Simulation
#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

def is_in_detector(hit_point: np.ndarray, detector_params: dict) -> bool:
    """Check if a photon is within a detector's volume."""
    detector_pos = detector_params["position"]
    detector_radius = detector_params["radius"]

    r = np.linalg.norm(hit_point[[0, 2]] - detector_pos[[0, 2]]) <= detector_radius
    y = detector_pos[1] <= hit_point[1] <= np.sign(detector_pos[1]) * np.abs(detector_pos[1] + detector_params["width"])

    return r and y

def detector_resolution(energy: float, detector_params: dict) -> float:
    return np.random.normal(energy, energy * detector_params["resolution"]) 

def exponential_bacground(energy: float, detector_params: dict) -> float:  
    """Simulate the exponential background."""
    return np.random.exponential(energy * detector_params["resolution"])

def gamma_detection(
        number_of_photons: int, 
        detector_params: dict, 
        step: float = 0.1, 
        detector_Z: int = Z, 
        source_position: np.ndarray = source_position
    ) -> float:

    """Simulate the detection of gamma rays by two detectors."""
    # photons_list = photons(number_of_photons)
    # photons_in_detector = []

    distance = np.linalg.norm(detector_params["position"] - source_position)
    
    # for photon in photons_list: 
    #     if photon.will_be_in_detector(detector_params, distance):
    #         photons_in_detector.append(photon)
    photons_in_detector = testing_photons(number_of_photons)
    photons_energies = []

    for photon in photons_in_detector:
        sigma_pe = cross_section_photoelectric(photon.energy)
        sigma_com = cross_section_compton(photon.energy)
        total_sigma = sigma_pe + sigma_com

        r1 = 0
        hit_point = photon.propagation(distance) 
        alpha_angle = np.arctan(hit_point[0] / hit_point[2]) if (hit_point[0] != 0 and hit_point[2] != 0) else 0
        h = distance - hit_point[1]  
        l = h / np.cos(alpha_angle)

        hit_point = photon.propagation(distance + l)

        while is_in_detector(hit_point, detector_params):
            if random.uniform(0, 1) < interaction_probability(photon, r1):
                if random.uniform(0, 1) < sigma_pe / total_sigma: # Photoelectric effect
                    photons_energies.append(detector_resolution(photon.energy, detector_params))
                else: # Compton scattering
                    # scattering_angle = random.uniform(0, np.pi)  # theta
                    scattering_angle = random_angles(photon)
                    phi = random.uniform(0, 2 * np.pi)  # phi
                    compton_energy = photon.energy - photon.compton_scattering(scattering_angle)
                    photons_energies.append(detector_resolution(compton_energy, detector_params))

                    # # Second Compton scattering
                    photon.energy = compton_energy
                    
                    x = np.sin(scattering_angle) * np.cos(phi)
                    y = np.sin(scattering_angle) * np.sin(phi)
                    z = np.cos(scattering_angle)
                    photon.direction = np.array([x, y, z])
                    r2 = 0
                    while is_in_detector(hit_point, detector_params): 
                        if random.uniform(0, 1) < interaction_probability(photon, r2):
                            if random.uniform(0, 1) < sigma_pe / total_sigma:
                                photons_energies.append(detector_resolution(photon.energy, detector_params))
                            else: 
                                scattering_angle = random.uniform(0, np.pi)  # theta
                                compton_energy = photon.energy - photon.compton_scattering(scattering_angle)
                                photons_energies.append(detector_resolution(compton_energy, detector_params))

                        r2 += step
                        hit_point += photon.propagation(r2)

            r1 += step
            hit_point = photon.propagation(distance + l + r1)
    return photons_energies                    

# Visualization
#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

def plot_energy_spectrum(energies, fileNamePNG, detector_name="Detector", bins=100, title="Energy Spectrum"):
    """Plot a histogram of the energy spectrum."""
    if not energies:
        print("No energies to plot.")
        return

    plt.figure(figsize=(10, 6))
    plt.hist(energies, bins=bins, color="blue", alpha=0.7, edgecolor="blue")
    
    plt.axvline(511, color='red', linestyle='--', linewidth=1.5, label='511 keV')
    plt.axvline(1274, color='red', linestyle='--', linewidth=1.5, label='1274 keV')
    plt.legend()

    plt.title(title, fontsize=16)
    plt.xlabel("Energy (keV)", fontsize=14)
    plt.ylabel("Counts", fontsize=14)
    plt.ylim(0, 250)
    plt.grid(alpha=0.4)
    plt.tight_layout()
    plt.savefig(fileNamePNG)
    plt.close()

# Example of usage
# number_of_photons = 10 ** 4 # Adjust the number of photons as needed
# energies_detected = gamma_detection(number_of_photons, detector_params["ugo"])
# plot_energy_spectrum(energies_detected, "ugo1.png", "Ugo", 1000, "")
# print("Fine\n")

# # Generate 10,000 random angles
# photon = Photon(511, np.array([0, 1, 0]))  # Initialize a photon object

# TEST PDF 
# # Generate random angles
# num_entries = 10000
# angles = [random_angles(photon) for _ in range(num_entries)]

# # Calculate PDF for comparison
# angles_pdf = np.linspace(0, np.pi, 500)
# pdf_values = [pdf_compton(photon, angle) for angle in angles_pdf]

# # Plot histogram of random angles
# plt.hist(angles, bins=50, density=True, alpha=0.5, label="Random Data (Histogram)")

# # Overlay the theoretical PDF
# plt.plot(angles_pdf, pdf_values, label="Compton PDF", color="blue", linewidth=2)

# # Labeling and formatting
# plt.xlabel("Angle (radians)")
# plt.ylabel("Probability Density")
# plt.title("Histogram of Random Angles vs. Compton PDF")
# plt.legend()
# plt.grid()
# plt.savefig("rn.png")

