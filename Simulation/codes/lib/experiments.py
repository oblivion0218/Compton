import numpy as np
import random
from . import particles as p
from . import detector as d
from . import source as s
from . import interactions as i
import matplotlib.patches as patches
from tqdm import tqdm

# Simulation pricipal axis is the y-axis

step = 0.1 #cm

def photon_propagation_to_target(photon: p.Photon, distance: float, direction=None) -> p.Photon:
    """
    Propagates a photon from the source to a target located at a specified distance and direction.

    :param photon: Photon object representing the gamma photon.
    :param distance: Distance between the source and the target (in cm).
    :param direction: Direction vector (x,y,z) pointing to the target. If None, y-axis [0,1,0] is used.
    :return: Photon object after propagation to the target.
    """
    # Set default direction to y-axis if None is provided
    if direction is None:
        direction = np.array([0, 1, 0])
    else:
        # Ensure direction is a numpy array and normalize it
        direction = np.array(direction)
        direction = direction / np.linalg.norm(direction)
    
    # Get the target position (a point on the target plane)
    target_position = photon.position + distance * direction
    
    # The normal vector of the target plane is the same as the direction to the target
    target_normal = -direction  # Negative because we want it facing toward the source
    
    # Calculate the vector from photon position to the target point
    vector_to_target = target_position - photon.position
    
    # Calculate the cosine of the angle between photon direction and target normal
    cos_angle = np.dot(photon.direction, target_normal)
    
    # Check if photon is moving toward the target
    if abs(cos_angle) < 1e-10:  # Nearly perpendicular, will never hit
        # Just propagate the original distance
        photon.propagation(distance)
        return photon
    
    # Calculate the distance to the intersection with the target plane
    # Using the plane equation: dot(target_normal, point - target_position) = 0
    propagation_distance = np.dot(target_normal, vector_to_target) / cos_angle
    
    # Propagate the photon to the intersection point
    if propagation_distance > 0:
        photon.propagation(propagation_distance)
    else:
        # If the intersection is behind the photon, just propagate the original distance
        photon.propagation(distance)
    
    return photon


def gamma_detection(photon: p.Photon, detector: d.Detector, distance_source_detector: float, step: float, true_energy = False) -> float:
    """
    Simulates the propagation and interaction of a gamma photon with a detector.
    
    :param photon: Photon object representing the gamma photon.
    :param detector: Detector object representing the detector.
    :param distance_source_detector: Distance between the source and the detector (in cm).
    :param step: Step size for photon propagation (in cm).
    :param true_energy: Flag to return the true energy of the electron or its detected energy.
    :return: Energy of the detected interaction in keV or 0 if no interaction occurs.
    """
    photon = photon_propagation_to_target(photon, distance_source_detector)

    # Initialize a variable for tracking the traveled distance within the detector
    electron = p.Electron(0, [0, 0, 0])
    # Propagate the photon within the detector until it exits or interacts
    while detector.is_inside(photon.position):
        # Check if the photon interacts with the detector material
        if random.uniform(0, 1) < i.interaction_probability(photon, step, detector):

            # Determine the interaction type (e.g., photoelectric or Compton)
            interaction = i.Interaction(i.which_interaction(photon, detector.Z))
            electron = interaction.interaction(photon)

            break  # Stop propagation after interaction

        photon.propagation(step)

    if true_energy == True: #If true_energy is True, returns the (real) electron's energy
        return electron.energy
    
    else: #If true_energy is False, returns the electron's energy after detection
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