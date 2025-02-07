import numpy as np
import particles as p
import matplotlib.patches as patches

# Constants
density = 3.67  # g/cm^3 (Density of NaI, the material of the detector)
Z_detector = 49.7  # Z_eff for NaI (effective atomic number of NaI material)
molar_mass = 149.89 # g/mol

class Detector:
    def __init__(
            self, 
            position: list[float], 
            radius: float, # cm
            width: float, #cm
            energetic_resolution: float, 
            Z: float = Z_detector, 
            density: float = density, # g/cm^3
            molar_mass: float = molar_mass # g/mol
        ):
        """
        Initialize the Detector object with parameters for position, size, and energetic resolution.

        :param position: List of 3 floats representing the detector's position in space [x, y, z].
        :param radius: The radius of the detector (in meters or appropriate units).
        :param width: The width of the detector (could represent depth or size depending on your setup).
        :param energetic_resolution: The energy resolution of the detector (could be in keV, % etc.).
        :param Z: atomic number of the detector material.
        :param density: density of the detector material.
        :param molar_mass: molar mass of the detector material.
        """
        self.position = np.array(position)  # Position of the detector in space (3D coordinates)
        self.radius = radius  # Radius of the detector (size)
        self.width = width  # Width of the detector (depth or size depending on setup)
        self.energetic_resolution = energetic_resolution  # Energy resolution of the detector (in keV or percentage)
        self.Z = Z  # Atomic number of the detector 
        self.density = density
        self.molar_mass = molar_mass 


    def __repr__(self):
        """
        String representation for the Detector object. This is called when printing the object.

        :return: A string representing the detector’s main properties.
        """
        return f"Detector(Position={self.position}, Radius={self.radius}, Width={self.width}, Energetic Resolution={self.energetic_resolution}, Z={self.Z}, Density={self.density}, Molar mass={self.molar_mass})"


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
        print(f"Density: {self.density}")
        print(f"Molar mass: {self.molar_mass}")

    
    def will_be_in_detector(self, point: np.ndarray, direction: np.ndarray) -> bool:
        """
        Determines if a particle will hit the detector based on its initial position and direction.

        :param point: The initial position of the particle (numpy array).
        :param direction: The direction vector of the particle (numpy array).
        :return: True if the particle will hit the detector, False otherwise.
        """
        # Calculate the distance between the particle's position and the detector's position
        distance = np.linalg.norm(point - self.position)
        
        # Calculate the angle of the detector's cone based on its radius and the distance
        theta_detector = np.arctan(self.radius / distance)
        
        # Calculate where the photon would hit the detector
        hit_point = point + distance * direction  # Intersection point of the particle's path and the detector
        
        # Calculate the distance in the X-Z plane between the hit point and the detector's position
        radius_hit = np.linalg.norm(hit_point[[0, 2]] - self.position[[0, 2]])
        
        # Calculate the angle at which the photon hits the detector
        theta_hit = np.arctan(radius_hit / hit_point[1])
        
        # Check if the photon is within the detector's radius and on the correct side
        within_radius = theta_hit <= theta_detector
        same_side = np.sign(hit_point[1]) == np.sign(self.position[1])  # Check if it's on the correct side of the detector
        
        return within_radius and same_side  # Return if both conditions are met
    

    def is_in_detector(self, point: np.ndarray) -> bool:
        """
        Checks if a given particle is inside the detector.

        :param point: The position of the particle (numpy array).
        :return: True if the particle is within the detector, False otherwise.
        """
        # Check if the particle is within the detector's radius in the X-Z plane
        r = np.linalg.norm(point[[0, 2]]) <= self.radius
        
        # Check if the Y-coordinate of the particle is within the bounds of the detector's width
        if self.position[1] > 0:
            y = self.position[1] <= point[1] <= self.position[1] + self.width
        else:
            y = self.position[1] >= point[1] >= self.position[1] - self.width

        return r and y  # Both conditions must be true to be inside the detector
    

    def resolution(self, energy: float) -> float:
        """
        Simulates the energy resolution of the detector by adding random fluctuations.

        :param energy: The true energy of the photon (float).
        :return: The energy after applying the detector's resolution (float).
        """
        return np.random.normal(energy, energy * self.energetic_resolution)  # Simulate energy resolution by adding Gaussian noise
    

    def detection(self, electron: p.Electron)-> float:
        """
        Simulates the energy detected by the interaction of an electron.

        :param electron: The incident electron.
        :return: The detected energy (simulated by applying resolution).
        """
        return self.resolution(electron.energy)  # Apply energy resolution to the detected electron's energy
    

    def draw_detector_3D(self, ax, axis='y', color='blue', alpha=0.5, label=None):
        """
        Draws the 3D representation of the detector using matplotlib.

        :param ax: The matplotlib axis object to plot on.
        :param axis: The axis along which the cylinder is oriented ('x', 'y', or 'z').
        :param color: Color of the cylinder.
        :param alpha: Transparency of the cylinder.
        :param label: Label for the detector in the legend.
        """
        # Create a set of theta values for plotting the circle
        theta = np.linspace(0, 2 * np.pi, 100)
        x = self.radius * np.cos(theta)
        y = self.radius * np.sin(theta)
    
        def linspace(base, width): 
            # Helper function to generate a linspace for the depth based on the position and width
            if base > 0:
                return np.linspace(base, base + width, 50)
            elif base < 0:
                return np.linspace(base, base - width, 50)
        
        # Create the grid for the cylinder's surface based on the specified axis
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
        
        # Plot the surface of the detector as a 3D object
        ax.plot_surface(X, Y, Z, color=color, alpha=alpha, edgecolor='none', label=label)
    
    
    def draw_detector_2D(self, ax, plane='xy', color='blue', label=None):
        """
        Draws the 2D projection of the detector on the specified plane.

        :param ax: The matplotlib axis object to plot on.
        :param plane: The plane on which to project the detector ('xy', 'xz', or 'yz').
        :param color: Color of the shape.
        :param label: Label for the shape in the legend.
        """
        if plane == 'xy':
            # Projection of the cylinder as a rectangle on the XY plane
            rect = patches.Rectangle(
                (self.position[0] - self.radius, self.position[1] if self.position[1] > 0 else self.position[1] - self.width),
                2 * self.radius,  # Rectangle width (diameter of the cylinder)
                self.width,      # Rectangle height (depth of the cylinder)
                color=color, alpha=0.5, label=label
            )
            ax.add_patch(rect)

        elif plane == 'xz':
            # Projection of the cylinder as a circle on the XZ plane
            circle = patches.Circle(
                (self.position[0], self.position[2]),
                self.radius,  # Cylinder radius
                color=color, alpha=0.5, label=label
            )
            ax.add_patch(circle)

        elif plane == 'yz':
            # Projection of the cylinder as a rectangle on the YZ plane
            rect = patches.Rectangle(
                (self.position[1] if self.position[1] > 0 else self.position[1] - self.width, self.position[2] - self.radius),
                self.width,  # Rectangle width (diameter of the cylinder)
                2 * self.radius,      # Rectangle height
                color=color, alpha=0.5, label=label
            )
            ax.add_patch(rect)

        else:
            raise ValueError("Plane must be 'xy', 'xz' or 'yz'.")

        # Adjust the axis limits to ensure proper visualization
        ax.set_aspect('equal', adjustable='datalim')
