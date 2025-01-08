import numpy as np
import random
import particles as p
import detector as d

# Physical constants
bn = 1e-28  # 1 barn in square meters, unit of cross-section
r_e = 2.817e-15  # Classical electron radius in meters
m_e = 511  # Electron rest mass energy in keV
alpha = 1 / 137  # Fine-structure constant (dimensionless)
 
# Cross section form "A Modern Primer to Particle and Nuclear Physycs" by F.Terranova page 69 - 71

# I don't consider the pair production because it's not relevant for the energy range we are considering
#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

# form "A Modern Primer to Particle and Nuclear Physycs" by F.Terranova page 69 - 71

def cross_section_photoelectric(photon: p.Photon, detector) -> float:
    """
    Calculate the photoelectric cross-section for a photon interacting with the detector material.

    :param photon: Photon object containing energy information.
    :param detector: Detector object containing material properties (e.g., atomic number Z).
    :return: Photoelectric cross-section in square meters.
    """
    epsilon = photon.energy / m_e  # Energy ratio between photon energy and electron rest mass energy
    if epsilon <= 1:  # Valid approximation for photon energy much less than electron rest mass energy
        c = 0.665 * np.sqrt(32) * alpha ** 4 * bn  # Empirical constant for low-energy range
        return c * detector.Z ** 5 / (epsilon ** 3.5)
    else:  # Valid approximation for photon energy much greater than electron rest mass energy
        c = 0.665 * (3 / 2) * alpha ** 4 * bn  # Empirical constant for high-energy range
        return c * detector.Z ** 5 / epsilon

def cross_section_compton(photon: p.Photon, detector) -> float:
    """
    Calculate the Compton scattering cross-section for a photon interacting with the detector material.

    :param photon: Photon object containing energy information.
    :param detector: Detector object containing material properties (e.g., atomic number Z).
    :return: Compton cross-section in square meters.
    """
    cross_section_thompson = (8 / 3) * np.pi * r_e ** 2  # Thomson scattering cross-section
    epsilon = photon.energy / m_e  # Energy ratio between photon energy and electron rest mass energy
    c = (3 / 4) * cross_section_thompson  # Scaling constant
    return c * detector.Z * (
        ((1 + epsilon) / epsilon ** 2) *
        ((2 * (1 + epsilon) / (1 + 2 * epsilon)) - (np.log(1 + 2 * epsilon) / epsilon)) +
        (np.log(1 + 2 * epsilon) / (2 * epsilon)) -
        (1 + 3 * epsilon) / (1 + 2 * epsilon) ** 2
    )

#-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
# 511: σ_pe = 2.4490838560315705e-28 ; σ_com = 1.3458389776243748e-27
# σ_pe / σ_tot = 0.15395806478424864 ; σ_com / σ_tot = 0.8460419352157514 (σ_tot = σ_pe + σ_com)
#-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
# 1274: σ_pe = 1.0008841411785894e-29 ; σ_com = 8.782021780069243e-28
# σ_pe / σ_tot = 0.01126854001241302 ; σ_com / σ_tot = 0.9887314599875869(σ_tot = σ_pe + σ_com)
#-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

# from wikipedia "https://en.wikipedia.org/wiki/Gamma_ray_cross_section"

#    Photoelectric: c = (16 / 3) * np.sqrt(2) * np.pi * r_e ** 2  * alpha ** 4
#    Compton: c = 2 * np.pi * r_e ** 2                                  

#-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
# 511: σ_pe = 1.2241772757589398e-28 ; σ_com = 1.345838977624375e-27
# σ_pe / σ_tot = 0.08337624282069687 ; σ_com / σ_tot = 0.9166237571793031 (σ_tot = σ_pe + σ_com)
#-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
# 1274: σ_pe = 5.002930456140891e-30 ; σ_com = 8.782021780069244e-28
# σ_pe / σ_tot = 0.005664517118619121 ; σ_com / σ_tot = 0.9943354828813808
#-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

# Interaction
#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

def interaction_probability(photon: p.Photon, distance: float) -> float:
    """
    Calculate the interaction probability for a photon over a given distance in the detector.

    :param photon: Photon object containing energy information.
    :param distance: Distance traveled in the detector (in cm).
    :return: Interaction probability as a dimensionless value.
    """
    attenuation_data = {  # Linear attenuation coefficients in cm^-1 for different photon energies
        10: 450,
        50: 25.0,
        100: 11.0,
        200: 4.5,
        511: 2.1,
        1000: 1.0,
        1274: 0.8
    }

    # Find the appropriate attenuation coefficient based on photon energy
    for attenuation_energy in attenuation_data:
        if photon.energy <= attenuation_energy:
            mu = attenuation_data[attenuation_energy]
            return 1 - np.exp(- mu * distance)

def which_interaction(photon: p.Photon, detector)-> str:
    sigma_pe = cross_section_photoelectric(photon, detector)
    sigma_com = cross_section_compton(photon, detector) 
    sigma_tot = sigma_pe + sigma_com

    r = random.uniform(0, 1)
    if r < sigma_pe / sigma_tot:
        return "photoelectric"
    else:
        return "compton"

def photoelectric_effect(photon: p.Photon)-> p.Electron:
    """
    Simulate the photoelectric effect, producing an electron with the same energy as the incident photon.

    :param photon: Incident photon object.
    :return: Electron object with energy equal to the incident photon energy.
    """
    energy = photon.energy
    photon.energy = 0
    return p.Electron(energy, photon.direction, photon.position)

def compton_scattering(photon: p.Photon, scattering_angle: float)-> p.Electron:
    """
    Simulate Compton scattering, producing a scattered photon and an electron.

    :param photon: Incident photon object.
    :return: Tuple containing the scattered photon and the electron.
    """
    theta = scattering_angle
    compton_energy_gamma = photon.compton_scattering(theta)
    compton_energy_electron = photon.energy - photon.compton_scattering(theta)  # Energy after scattering

    photon.energy = compton_energy_gamma  # Update photon energy

    phi = random.uniform(0, 2 * np.pi)
    electron_direction = np.array([np.sin(scattering_angle) * np.cos(phi), np.sin(scattering_angle) * np.sin(phi),  np.cos(scattering_angle)])

    return p.Electron(compton_energy_electron, electron_direction, photon.position)  # Create an electron with the transferred energy

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
        print("interaction:" + self.type)
    
    def cross_section(self, photon: p.Photon, detector) -> float:
        """
        Calculate the cross-section for the specified interaction type.

        :param photon: Photon object containing energy information.
        :param detector: Detector object containing material properties (e.g., atomic number Z).
        :return: Cross-section for the interaction in square meters.
        """
        if self.type == "photoelectric":
            return cross_section_photoelectric(photon, detector)
        if self.type == "compton":
            return cross_section_compton(photon, detector)

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
