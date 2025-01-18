import numpy as np
import random
import particles as p

# Physical constants
r_e = 2.817e-13  # Classical electron radius in cm
m_e = 511  # Electron rest mass energy in keV
alpha = 1 / 137  # Fine-structure constant (dimensionless)
N_a = 6.022e23 # mol^(-1)
 
# Cross section from wikipedia "https://en.wikipedia.org/wiki/Gamma_ray_cross_section"

# I don't consider the pair production because it's not relevant for the energy range we are considering
#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

def cross_section_thomson():
    """
    Calculate the Thomson cross-section, which represents the scattering 
    of electromagnetic waves by free electrons at low photon energies.

    :return: Thomson cross-section in square meters.
    """
    return (8 / 3) * np.pi * r_e ** 2


def cross_section_photoelectric(photon: p.Photon, Z: float) -> float:
    """
    Calculate the photoelectric cross-section for a photon interacting with a material.

    :param photon: Photon object containing energy information (energy in MeV).
    :param Z: Atomic number of the material.
    :return: Photoelectric cross-section in square meters.
    """
    bond_energy = 0.007  # Example binding energy in MeV
    gamma = (photon.energy + m_e - bond_energy) / m_e  # Lorentz factor
    c = (3 / 2) * (alpha ** 4) * cross_section_thomson()  # Coefficient based on fine-structure constant

    return c * (
        ((Z * m_e / photon.energy) ** 5) *
        (gamma ** 2 - 1) ** (3/2) *
        (
            (4/3) + (gamma * (gamma - 2) / (gamma + 1)) *
            (
                1 - (1 / (2 * gamma * (gamma ** 2 - 1)**(1/2))) *
                np.log((gamma + (gamma ** 2 - 1)**(1/2)) / (gamma - (gamma ** 2 - 1)**(1/2)))
            )
        )
    )


def cross_section_compton(photon: p.Photon, Z: float) -> float:
    """
    Calculate the Compton scattering cross-section for a photon interacting with a material.

    :param photon: Photon object containing energy information (energy in MeV).
    :param Z: Atomic number of the material.
    :return: Compton cross-section in square meters.
    """
    epsilon = photon.energy / m_e  # Ratio of photon energy to electron rest mass energy
    c = (3 / 4) * cross_section_thomson()  # Coefficient based on Thomson cross-section

    return c * Z * (
        ((1 + epsilon) / epsilon ** 2) *
        ((2 * (1 + epsilon) / (1 + 2 * epsilon)) - (np.log(1 + 2 * epsilon) / epsilon)) +
        (np.log(1 + 2 * epsilon) / (2 * epsilon)) -
        ((1 + 3 * epsilon) / (1 + 2 * epsilon) ** 2)
    )

# For Z = 49.7

# Energy: 511 keV
# σ_pe = 5.5397e-28 ; σ_com = 1.4031e-27
# σ_pe / σ_tot = 0.2831 ; σ_com / σ_tot = 0.7169 (σ_tot = 1.9571e-27)
# --------------------------------------------------------------------------------
# Energy: 1274 keV
# σ_pe = 7.1684e-29 ; σ_com = 9.1557e-28
# σ_pe / σ_tot = 0.0726 ; σ_com / σ_tot = 0.9274 (σ_tot = 9.8726e-28)
# --------------------------------------------------------------------------------


# Interaction
#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

def attenuation_factor(photon: p.Photon, distance: float, total_cross_section: float, scattering_target) -> float:
    """
    Calculate the attenuation factor for a photon traveling through a material.

    :param photon: Photon object containing energy information.
    :param distance: Distance traveled in the material (in cm).
    :param total_cross_section: Total cross-section for the interaction (in cm²).
    :param scattering_target: Material properties (e.g., density, molar mass).
    :return: Attenuation factor (inverse mean free path).
    """
    number_of_scattering_centers = (scattering_target.density * N_a) / scattering_target.molar_mass
    return total_cross_section * number_of_scattering_centers # cm^-1


def interaction_probability(photon: p.Photon, distance: float, scattering_target) -> float:
    """
    Calculate the probability of interaction for a photon traveling a given distance in a material.

    :param photon: Photon object containing energy information.
    :param distance: Distance traveled in the material (in cm).
    :param total_cross_section: Total cross-section for the interaction (in cm²).
    :param scattering_target: Material properties (e.g., density, molar mass).
    :return: Interaction probability as a dimensionless value.
    """
    total_cross_section = cross_section_photoelectric(photon, scattering_target.Z) + cross_section_compton(photon, scattering_target.Z)

    mu = attenuation_factor(photon, distance, total_cross_section, scattering_target)

    return 1 - np.exp(-mu * distance)


def which_interaction(photon: p.Photon, Z: float) -> str:
    """
    Determine the type of interaction (photoelectric or Compton) for a photon.

    :param photon: Photon object containing energy information.
    :param Z: Atomic number of the material.
    :return: A string indicating the interaction type ("photoelectric" or "compton").
    """
    sigma_pe = cross_section_photoelectric(photon, Z)
    sigma_com = cross_section_compton(photon, Z)
    sigma_tot = sigma_pe + sigma_com

    r = random.uniform(0, 1)
    return "photoelectric" if r < sigma_pe / sigma_tot else "compton"


def photoelectric_effect(photon: p.Photon) -> p.Electron:
    """
    Simulate the photoelectric effect, where a photon transfers all its energy to an electron.

    :param photon: Incident photon object.
    :return: Electron object with energy equal to the incident photon's energy.
    """
    energy = photon.energy
    photon.energy = 0  # Photon is fully absorbed
    return p.Electron(energy, photon.direction, photon.position)


def compton_scattering(photon: p.Photon, scattering_angle: float) -> p.Electron:
    """
    Simulate Compton scattering, producing a scattered photon and an electron.

    :param photon: Incident photon object.
    :param scattering_angle: Angle of photon scattering (in radians).
    :return: Scattered electron object with energy from the interaction.
    """
    theta = scattering_angle
    compton_energy_gamma = photon.compton_scattering(theta)
    compton_energy_electron = photon.energy - compton_energy_gamma  # Energy transferred to the electron

    photon.energy = compton_energy_gamma  # Update photon energy

    phi = random.uniform(0, 2 * np.pi)  # Azimuthal angle
    electron_direction = np.array([
        np.sin(scattering_angle) * np.cos(phi),
        np.sin(scattering_angle) * np.sin(phi),
        np.cos(scattering_angle)
    ])

    return p.Electron(compton_energy_electron, electron_direction, photon.position)


class Interaction:
    def __init__(self, type: str):
        """
        Initialize an interaction type (e.g., photoelectric or Compton).

        :param type: Type of interaction as a string ("photoelectric" or "compton").
        """
        self.type = type

    
    def __repr__(self):
        """
        Return a string representation of the interaction type.

        :return: Interaction type as a string.
        """
        return self.type 


    def info(self):
        """
        Print information about the interaction.

        :return: None.
        """
        print(f"Interaction: {self.type}")


    def cross_section(self, photon: p.Photon, Z) -> float:
        """
        Calculate the cross-section for the specified interaction type.

        :param photon: Photon object containing energy information.
        :param detector: Detector object containing material properties (e.g., atomic number Z).
        :return: Cross-section for the interaction in square meters.
        """
        if self.type == "photoelectric":
            return cross_section_photoelectric(photon, Z)
        if self.type == "compton":
            return cross_section_compton(photon, Z)


    def interaction(self, photon: p.Photon)-> p.Electron:
        """
        Perform the specified interaction, returning the resulting products.

        :param photon: Photon object containing energy information.
        :return: Result of the interaction (electron for photoelectric or photon-electron pair for Compton).
        """
        if self.type == "photoelectric":
            return photoelectric_effect(photon)
        if self.type == "compton":
            return compton_scattering(photon, photon.compton_angle())