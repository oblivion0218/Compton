import numpy as np
import random
import particles as p
import detector as d
import source as s
import interactions as i
import matplotlib.patches as patches
        
step = 0.1 #cm

def photon_propagation_to_target(photon: p.Photon, distance_source_detector: float) -> p.Photon: # figure: geometry_exp.png
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
    while detector.is_in_detector(photon.position):

        total_cross_section = i.cross_section_photoelectric(photon, detector.Z) + i.cross_section_compton(photon, detector.Z)

        # Check if the photon interacts with the detector material
        if random.uniform(0, 1) < i.interaction_probability(photon, step, detector):

            # Determine the interaction type (e.g., photoelectric or Compton)
            interaction = i.Interaction(i.which_interaction(photon, detector.Z))
            electron = interaction.interaction(photon)

            break  # Stop propagation after interaction

        photon.propagation(step)

    return detector.detection(electron)


def spectroscopy_measurement(number_of_photons, detector: d.Detector, testing: bool = False, step: float = step) -> list[float]:
    """
    Simulates the interaction of multiple gamma photons with a detector to calculate detected energies.
    
    :param number_of_photons: Number of photons to simulate.
    :param detector: Detector object where photons are detected.
    :param testing: Flag to enable testing mode, which uses predefined photons.
    :param step: Step size for photon propagation (in cm).
    :return: List of detected photon energies (in keV).
    """
    source = s.Source()
    direction = [0, np.sign(detector.position[1]), 0]
    # Generate photons either for testing or normal emission
    photons = source.testing_photons(number_of_photons, direction) if testing else source.photon_emission(number_of_photons)

    distance = np.linalg.norm(detector.position - source.position)
    detected_energies = []

    for photon in photons: 
        energy = gamma_detection(photon, detector, distance, step)
        detected_energies.append(energy)
    
    return detected_energies


def coincidence_photons(number_of_photons: int, gate_detector: d.Detector, spettroscopy_detector: d.Detector, testing: bool = False, step: float = step)-> list[p.Photon]:
    """
    Simulates the interaction of multiple gamma photons with two detectors to calculate detected energies.
    
    :param number_of_photons: Number of photons to simulate.
    :param gate_detector: Detector object for the gate detector.
    :param spettroscopy_detector: Detector object for the spectroscopy detector.
    :param testing: Flag to enable testing mode, which uses predefined photons.
    :param step: Step size for photon propagation (in cm).
    :return: List of detected photon energies (in keV).
    """
    energies = {511: 1, 1274: 0}
    source = s.Source(energies)
    direction = [0, np.sign(gate_detector.position[1]), 0]
    # Generate photons either for testing or normal emission
    photons = source.testing_photons(number_of_photons, direction) if testing else source.photon_emission(number_of_photons)

    distance_gate = np.linalg.norm(gate_detector.position - source.position)
    distance_spettroscopy = np.linalg.norm(spettroscopy_detector.position - source.position)
    coincidence_photons = []

    for photon in photons:
        photon_energy = photon.energy
        energy_gate = gamma_detection(photon, gate_detector, distance_gate, step)

        if energy_gate > 0:
            spettroscopy_photon = p.Photon(photon_energy, (-1) * photon.direction)
            coincidence_photons.append(spettroscopy_photon)
    
    return coincidence_photons


def coincidence_measurement(number_of_photons: int, gate_detector: d.Detector, spettroscopy_detector: d.Detector, testing: bool = False, step: float = step) -> list[float]:
    """
    Simulates the detection of coincident gamma photons in two detectors.
    
    :param number_of_photons: Number of photons to simulate.
    :param gate_detector: Detector object for the gate detector.
    :param spettroscopy_detector: Detector object for the spectroscopy detector.
    :param testing: Flag to enable testing mode, which uses predefined photons.
    :param step: Step size for photon propagation (in cm).
    :return: List of detected photon energies (in keV).
    """
    photons = coincidence_photons(number_of_photons, gate_detector, spettroscopy_detector, testing, step)
    detected_energies = []
    distance_spettroscopy = np.linalg.norm(spettroscopy_detector.position) # source in [0, 0, 0]

    for photon in photons:
        energy_spettroscopy = gamma_detection(photon, spettroscopy_detector, distance_spettroscopy, step)
        detected_energies.append(energy_spettroscopy)

    return detected_energies
