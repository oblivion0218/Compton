import numpy as np
import particles as p

density = 3.67  # g/cm^3 (Density of NaI, detector material)
Z_detector = 47  # Approximate atomic number of NaI (mean of Na and I)

class Detector:
    def __init__(self, position: list[float], radius: float, width: float, energetic_resolution: float, Z: int = Z_detector):
        """
        Initialize the Detector object with parameters for position, size, and energetic resolution.

        :param position: List of 3 floats representing the detector's position in space [x, y, z].
        :param radius: The radius of the detector (in meters or appropriate units).
        :param width: The width of the detector (this could represent depth or size depending on your setup).
        :param energetic_resolution: The energy resolution of the detector (this could be in keV, % etc.)
        :param Z: atomic number of the detector material
        """
        self.position = np.array(position)  # Position of the detector in space
        self.radius = radius  # Detector's radius (size)
        self.width = width  # Detector's width (could represent depth)
        self.energetic_resolution = energetic_resolution  # Detector's energy resolution (in keV or percentage)
        self.Z = Z

    def __repr__(self):
        """
        String representation for the Detector object. This is called when printing the object.

        :return: A string representing the detector’s main properties.
        """
        return f"Detector(Position={self.position}, Radius={self.radius}, Width={self.width}, Energetic Resolution={self.energetic_resolution}, Z={self.Z})"

    def info(self):
        """
        Print the detailed information of the detector.

        :return: None
        """
        print("Detector Information:")
        print(f"Position: {self.position}")
        print(f"Radius: {self.radius} m")
        print(f"Width: {self.width} m")
        print(f"Energetic Resolution: {self.energetic_resolution} keV")
        print(f"Z: {self.Z}")
    
    def will_be_in_detector(self, point: np.ndarray, direction: np.ndarray) -> bool:
        """
        Determines if a particle will hit the detector based on its initial position and direction.

        :param point: The initial position of the particle (numpy array).
        :param direction: The direction vector of the particle (numpy array).
        :return: True if the particle will hit the detector, False otherwise.
        """
        distance = np.linalg.norm(point - self.position)  # Calculate the distance from the point to the detector
        theta_detector = np.arctan(self.radius / distance)  # Angle of the detector cone         
        
        hit_point = point + distance * direction  # Calculate where the photon would hit
        radius_hit = np.linalg.norm(hit_point[[0, 2]] - self.position[[0, 2]])  # Distance in the x-z plane
        theta_hit = np.arctan(radius_hit / hit_point[1])  # Angle at which the photon hits         
        
        # Check if the photon is within the detector's radius and on the correct side
        within_radius = theta_hit <= theta_detector
        same_side = np.sign(hit_point[1]) == np.sign(self.position[1])  # Check if the photon is on the correct side         
        
        return within_radius and same_side  # Return if both conditions are met
    
    def is_in_detector(self, point: np.ndarray) -> bool:
        """
        Checks if a given particle is inside the detector.

        :param point: The position of the particle (numpy array).
        :return: True if the particle is within the detector, False otherwise.
        """
        r = np.linalg.norm(point[[0, 2]]) <= self.radius  # Check if the photon is within the detector's radius in x-z plane
        y = self.position[1] <= point[1] <= np.sign(self.position[1]) * np.abs(self.position[1] + self.width)  # Check if y is within the bounds         
        
        return r and y  # Both conditions must be true
    
    def resolution(self, energy: float) -> float:
        """
        Simulates the energy resolution of the detector by adding random fluctuations.

        :param energy: The true energy of the photon (float).
        :return: The energy after applying the detector's resolution (float).
        """
        return np.random.normal(energy, energy * self.energetic_resolution)  # Simulate energy resolution         
    
    def detection(self, electron: p.Electron)-> float:
        """
        Simulates the energy detected by the interaction of an electron.

        :param electron: The incident electron.
        :return: The detected energy.
        """
        return self.resolution(electron.energy)