import numpy as np
import random
from . import particles as p
from . import detector as d
from . import source as s
from . import interactions as i
import matplotlib.patches as patches
from tqdm import tqdm

step = 0.1 #cm

def photon_propagation_to_target(photon: p.Photon, distance_source_detector: float) -> p.Photon: # figure: geometry_exp.png
    """
    Propagates a photon from the source to the detector.

    :param photon: Photon object representing the gamma photon.
    :param distance_source_detector: Distance between the source and the detector (in cm).
    :return: Photon object after propagation to the detector.
    """
    photon.propagation(distance_source_detector)
    # Calculate the angle (α) between the photon's trajectory and the detector surface
    alpha_angle = np.arctan(photon.position[0] / photon.position[2]) if (photon.position[0] != 0 and photon.position[2] != 0) else 0
    
    # Calculate the vertical offset (H) and the diagonal length (L) to the detector
    H = distance_source_detector - np.sign(photon.position[1]) * photon.position[1]
    L = H / np.cos(alpha_angle)

    # Adjust the propagation distance based on the offset
    photon.propagation(L)

    return photon


def gamma_detection(photon: p.Photon, detector: d.Detector, distance_source_detector: float, step: float) -> float:
    """
    Simulates the propagation and interaction of a gamma photon with a detector.
    
    :param photon: Photon object representing the gamma photon.
    :param detector: Detector object representing the detector.
    :param distance_source_detector: Distance between the source and the detector (in cm).
    :param step: Step size for photon propagation (in cm).
    :return: Energy of the detected interaction in keV or 0 if no interaction occurs.
    """
    photon = photon_propagation_to_target(photon, distance_source_detector)

    # Initialize a variable for tracking the traveled distance within the detector
    electron = p.Electron(0, [0, 0, 0])
    # Propagate the photon within the detector until it exits or interacts
    while detector.is_inside(photon.position):
        total_cross_section = i.cross_section_photoelectric(photon, detector.Z) + i.cross_section_compton(photon, detector.Z)

        # Check if the photon interacts with the detector material
        if random.uniform(0, 1) < i.interaction_probability(photon, step, detector):

            # Determine the interaction type (e.g., photoelectric or Compton)
            interaction = i.Interaction(i.which_interaction(photon, detector.Z))
            electron = interaction.interaction(photon)

            break  # Stop propagation after interaction

        photon.propagation(step)

    return detector.detection(electron)


def spectroscopy_measurement(number_of_photons, detector: d.Detector, source: s.Source, testing: bool = False, step: float = step) -> list[float]:
    """
    Simulates the interaction of multiple gamma photons with a detector to calculate detected energies.
    
    :param number_of_photons: Number of photons to simulate.
    :param detector: Detector object where photons are detected.
    :param testing: Flag to enable testing mode, which uses predefined photons.
    :param step: Step size for photon propagation (in cm).
    :return: List of detected photon energies (in keV).
    """
    center_detector = detector.center()
    len_principal_axis = np.linalg.norm(detector.principal_axis())
    direction = [0, np.sign(center_detector[1]), 0]
    # Generate photons either for testing or normal emission
    photons = source.testing_photons(number_of_photons, direction) if testing else source.photon_emission(number_of_photons)
    distance = np.linalg.norm(center_detector - source.position) - len_principal_axis/2
    detected_energies = []

    for photon in tqdm(photons, desc="Simulating photons", unit="photon"): 
        energy = gamma_detection(photon, detector, distance, step)
        detected_energies.append(energy)
    
    return detected_energies