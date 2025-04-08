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

def photon_propagation_to_target(photon: p.Photon, target: d.Object) -> p.Photon:
    """
    Propagates a photon from its current position to the first face of a target object.
    
    :param photon: Photon object representing the gamma photon.
    :param target: Target object (instance of d.Object) that the photon will propagate toward.
    :return: Photon object after propagation to the target's first face.
    """
    # Get the first face position (position[0])
    first_face_position = np.array(target.position[0])
    
    # Calculate principal axis of the object
    principal_axis = target.principal_axis()
    
    # The normal vector of the first face is opposite to the principal axis (facing toward the source)
    target_normal = -principal_axis / np.linalg.norm(principal_axis)
    
    # Calculate the vector from photon position to the first face
    vector_to_face = first_face_position - photon.position
    
    # Calculate the cosine of the angle between photon direction and target normal
    cos_angle = np.dot(photon.direction, target_normal)
    
    # Check if photon is moving toward the target face
    if abs(cos_angle) < 1e-10:  # Nearly perpendicular, will never hit
        # Just propagate toward the first face position
        distance_to_face = np.linalg.norm(vector_to_face)
        photon.propagation(distance_to_face)
        return photon
    
    # Calculate the distance to the intersection with the target plane
    # Using the plane equation: dot(target_normal, point - target_position) = 0
    propagation_distance = np.dot(target_normal, vector_to_face) / cos_angle
    
    # Propagate the photon to the intersection point
    if propagation_distance > 0:
        photon.propagation(propagation_distance)
        
        # Check if the intersection point is within the radius of the first face
        intersection_point = photon.position
        
        # Calculate the vector from the first face center to the intersection point
        vector_on_face = intersection_point - first_face_position
        
        # Project this vector onto the plane
        projection = vector_on_face - np.dot(vector_on_face, target_normal) * target_normal
        
        # Check if the projection length is less than the radius
        if np.linalg.norm(projection) <= target.radius:
            return photon  # The photon hit the face within its radius
        else:
            # If missed the circular face, propagate additional distance
            distance_to_face = np.linalg.norm(first_face_position - photon.position)
            photon.propagation(distance_to_face)
    else:
        # If the intersection is behind the photon, just propagate toward the face
        distance_to_face = np.linalg.norm(vector_to_face)
        photon.propagation(distance_to_face)
    
    return photon


def gamma_detection(photon: p.Photon, detector: d.Detector, step: float, true_energy = False) -> float:
    """
    Simulates the propagation and interaction of a gamma photon with a detector.
    
    :param photon: Photon object representing the gamma photon.
    :param detector: Detector object representing the detector.
    :param step: Step size for photon propagation (in cm).
    :param true_energy: Flag to return the true energy of the electron or its detected energy.
    :return: Energy of the detected interaction in keV or 0 if no interaction occurs.
    """
    photon = photon_propagation_to_target(photon, detector)

    if photon.energy <= 0:
        return 0

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
    detected_energies = []

    for photon in tqdm(photons, desc="Simulating photons", unit="photon"): 
        energy = gamma_detection(photon, detector, step)
        detected_energies.append(energy)
    
    return detected_energies