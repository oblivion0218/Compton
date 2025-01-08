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
    distance_source_detector += L
    photon.propagation(distance_source_detector)

    return photon

# Function to simulate the detection of a single gamma photon
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
    r = 0
    electron = p.Electron(0, [0, 0, 0])
    # Propagate the photon within the detector until it exits or interacts
    while detector.is_in_detector(photon.position):

        # Check if the photon interacts with the detector material
        if random.uniform(0, 1) < i.interaction_probability(photon, r):

            # Determine the interaction type (e.g., photoelectric or Compton)
            interaction = i.Interaction(i.which_interaction(photon, detector))
            electron = interaction.interaction(photon)

            break  # Stop propagation after interaction

        r += step
        photon.propagation(distance_source_detector + r)


    return detector.detection(electron)

# Function to simulate the detection of multiple gamma photons
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

Z_Pb = 82
density_Pb = 11.34

class Target:
    def __init__(self, position: list[float], radius: float, width: float, density: float = 11.34, Z: float = 82):
        """
        Initialize the target with its position, dimensions, density, and atomic number.
        
        :param position: Position of the target as a 3D vector (list of floats).
        :param radius: Radius of the target in meters.
        :param width: Width (height) of the target in meters.
        :param density: Density of the target material in g/cm^3 (default: density of Pb).
        :param Z: Atomic number of the target material (default: Z of Pb).
        """
        self.position = np.array(position)
        self.radius = radius
        self.width = width
        self.density = density
        self.Z = Z

    def __repr__(self):
        """
        String representation for the Target object.
        """
        return f"Target(Position={self.position}, Radius={self.radius}, Width={self.width}, Density={self.density}, Z={self.Z})"
    
    def info(self):
        """
        Print the detailed information of the target.
        """
        print("Target Information:")
        print(f"Position: {self.position}")
        print(f"Radius: {self.radius} m")
        print(f"Width: {self.width} m")
        print(f"Density: {self.density} g/cm^3")
        print(f"Z: {self.Z}")
    
    def is_in_target(self, point: np.ndarray) -> bool:
        """
        Checks if a given point (position of a photon) is inside the target.
        
        :param point: The position of the particle (numpy array).
        :return: True if the point is within the target, False otherwise.
        """
        r = np.linalg.norm(point[[0, 2]]) <= self.radius  # Check radial distance
        if self.position[1] > 0:
            y = self.position[1] <= point[1] <= self.position[1] + self.width  # Check y-range
        else:
            y = self.position[1] >= point[1] >= self.position[1] - self.width
        
        return r and y
    
    def attenuation(self, photon: p.Photon, distance: float) -> float:
        """
        Calculate the attenuation of a gamma photon passing through the target.
        
        :param photon: Photon object with energy and position.
        :param distance: Distance the photon travels through the material.
        :return: Attenuation factor.
        """
        # Linear attenuation coefficients for different photon energies (in cm^-1)
        attenuation_data = {100: 65, 511: 3.2, 1274: 0.7}

        # Find the appropriate attenuation coefficient based on photon energy
        for attenuation_energy in attenuation_data:
            if photon.energy <= attenuation_energy:
                mu = attenuation_data[attenuation_energy]
                return 1 - np.exp(-mu * distance)  # Attenuation formula

    def draw_target_3D(self, ax, axis='y', color='grey', alpha=0.5, label=None):
        """
        Draw the target as a cylinder in 3D space, projected along a specified axis.
        
        :param ax: The 3D axis where the target will be drawn.
        :param axis: The axis ('x', 'y', or 'z') along which the cylinder is oriented.
        :param color: Color of the target.
        :param alpha: Transparency of the target.
        :param label: Label for the target (optional).
        """
        theta = np.linspace(0, 2 * np.pi, 100)
        x = self.radius * np.cos(theta)
        y = self.radius * np.sin(theta)
    
        def linspace(base, width): 
            """Generate linear space for positioning in 3D."""
            return np.linspace(base, base + width, 50) if base > 0 else np.linspace(base, base - width, 50)
        
        if axis == 'z':
            z = linspace(self.position[2], self.width)
            X, Z = np.meshgrid(x, z)
            Y, _ = np.meshgrid(y, z)
            X += self.position[0]
            Y += self.position[1]
        elif axis == 'x':
            z = linspace(self.position[0], self.width)
            Z, X = np.meshgrid(x, z)
            Y, _ = np.meshgrid(y, z)
            Z += self.position[2]
            Y += self.position[1]
        elif axis == 'y':
            z = linspace(self.position[1], self.width)
            Z, Y = np.meshgrid(x, z)
            X, _ = np.meshgrid(y, z)
            Z += self.position[2]
            X += self.position[0]
        else:
            raise ValueError("Axis must be 'x', 'y' or 'z'")
        
        # Color the surface of the cylinder
        ax.plot_surface(X, Y, Z, color=color, alpha=alpha, edgecolor='none', label=label)

    def draw_target_2D(self, ax, plane='xy', color='grey', label=None):
        """
        Draw the 2D projection of the target on the specified plane.
        
        :param ax: The axis where the 2D projection will be drawn.
        :param plane: The plane ('xy', 'xz', or 'yz') for the projection.
        :param color: Color of the projection.
        :param label: Label for the projection (optional).
        """
        if plane == 'xy':
            rect = patches.Rectangle(
                (self.position[0] - self.radius, self.position[1] if self.position[1] > 0 else self.position[1] - self.width),
                2 * self.radius,  # Width of the rectangle (diameter of the cylinder)
                self.width,       # Height of the rectangle (depth of the cylinder)
                color=color, alpha=0.5, label=label
            )
            ax.add_patch(rect)
        elif plane == 'xz':
            circle = patches.Circle(
                (self.position[0], self.position[2]),
                self.radius,  # Radius of the cylinder
                color=color, alpha=0.5, label=label
            )
            ax.add_patch(circle)
        elif plane == 'yz':
            rect = patches.Rectangle(
                (self.position[1] if self.position[1] > 0 else self.position[1] - self.width, self.position[2] - self.radius),
                self.width,  # Width of the rectangle (diameter of the cylinder)
                2 * self.radius,  # Height of the rectangle
                color=color, alpha=0.5, label=label
            )
            ax.add_patch(rect)
        else:
            raise ValueError("Plane must be 'xy', 'xz' or 'yz'.")
        
        ax.set_aspect('equal', adjustable='datalim')

def target_scattering(number_of_photons: int, scattering_angle: float, target: Target, gate_detector: d.Detector, spettroscopy_detector: d.Detector, testing: bool = False, step: float = 0.1) -> list[p.Photon]:
    """
    Simulates the scattering of photons in a target and detection by detectors.
    
    :param number_of_photons: Number of photons to simulate.
    :param scattering_angle: The scattering angle (in radians).
    :param target: The target material object.
    :param gate_detector: Gate detector object.
    :param spettroscopy_detector: Spectroscopy detector object.
    :param testing: If True, use predefined photons for testing.
    :param step: Step size for photon propagation (in cm).
    :return: List of photons detected after scattering.
    """
    energies = {511: 1, 1274: 0}
    source = s.Source(energies)
    direction = [0, np.sign(gate_detector.position[1]), 0]
    # Generate photons either for testing or normal emission
    photons = source.testing_photons(number_of_photons, direction) if testing else source.photon_emission(number_of_photons)

    distance_gate_source = np.linalg.norm(gate_detector.position - source.position)
    distance_source_target = np.linalg.norm(target.position - source.position)
    scattering_photons = []

    for photon in photons:
        photon_energy = photon.energy
        energy_gate = gamma_detection(photon, gate_detector, distance_gate_source, step)
        
        if energy_gate > 0:
            sibling_photon = p.Photon(photon_energy, (-1) * photon.direction)
            sibling_photon = photon_propagation_to_target(sibling_photon, distance_source_target)

            r = 0
            scattering_photon = p.Photon(0, [0, 0, 0])
            while target.is_in_target(sibling_photon.position):
                if random.uniform(0, 1) < i.interaction_probability(sibling_photon, r):
                    interaction = i.Interaction(i.which_interaction(sibling_photon, target))
                    if interaction.type == "compton":
                        scattering_photon.energy = sibling_photon.compton_scattering(scattering_angle)
                        phi = random.uniform(0, 2 * np.pi)

                        scattering_photon.direction = np.array([np.sin(scattering_angle) * np.sin(phi), np.cos(scattering_angle), np.sin(scattering_angle) * np.cos(phi)])
                        
                        scattering_photon.position = sibling_photon.position

                        scattering_photons.append(scattering_photon)
                    break

                r += step
                sibling_photon.propagation(distance_source_target + r)

    # Manca Photon attenuation check
    return scattering_photons

def compton_measurement(number_of_photons: int, scattering_angle: float, target: Target, gate_detector: d.Detector, spettroscopy_detector: d.Detector, testing: bool = False, step: float = 0.1) -> list[float]:
    """
    Measures the energy of photons detected after Compton scattering.
    
    :param number_of_photons: Number of photons to simulate.
    :param scattering_angle: The scattering angle in radians.
    :param target: Target object.
    :param gate_detector: Gate detector object.
    :param spettroscopy_detector: Spectroscopy detector object.
    :param testing: If True, uses predefined photons for testing.
    :param step: Step size for photon propagation.
    :return: List of detected energies.
    """
    scattering_photons = target_scattering(number_of_photons, scattering_angle, target, gate_detector, spettroscopy_detector, testing, step)
    detected_energies = []
    distance_spettroscopy_target = np.linalg.norm(spettroscopy_detector.position - target.position)

    for photon in scattering_photons:
        energy_spettroscopy = gamma_detection(photon, spettroscopy_detector, distance_spettroscopy_target, step)
        detected_energies.append(energy_spettroscopy)

    return detected_energies
