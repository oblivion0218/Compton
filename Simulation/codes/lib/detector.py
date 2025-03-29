import numpy as np
import particles as p
import matplotlib.patches as patches

class Object:
    def __init__(self, position: tuple[list[float], list[float]], radius: float, Z: float, density: float, molar_mass: float):        
        """
        Initialize the object with its position and dimensions.

        :param position: Position of the object as a 3D vector (list of floats), the center of the two faces.
        :param radius: Radius of the object in meters.
        :param Z: Atomic number of the object.
        :param density: Density of the object material in g/cm^3.
        :param molar_mass: Molar mass of the object material in g/mol.
        """
        self.position = position
        self.radius = radius
        self.Z = Z  # Atomic number of the object
        self.density = density
        self.molar_mass = molar_mass

    def __repr__(self):
        """
        String representation for the Object object.

        :return: A string representing the object's main properties.
        """
        return f"Object(Position={self.position}, Radius={self.radius}, Z={self.Z}, Density={self.density}, Molar mass={self.molar_mass})"
    

    def info(self):
        """
        Print the detailed information of the object.

        :return: None
        """
        print("Object Information:")
        print(f"Position: {self.position}")
        print(f"Radius: {self.radius} m")
        print(f"Z: {self.Z}")
        print(f"Density: {self.density}")
        print(f"Molar mass: {self.molar_mass}")


    def is_inside(self, point: np.ndarray) -> bool:
        """
        Checks if a given point (position of a photon) is inside the object.

        :param point: The position of the particle (numpy array).
        :return: True if the point is within the object, False otherwise.
        """
        # Define the principal axis of the object, the one that is perpendicular to the two faces
        # This axis is the one that connect the two point inside the position tuple
        principal_axis = np.array(self.position[1]) - np.array(self.position[0])

        # Rotate the axis (x, y, z) to align z with the principal axis, and the origin with the middle of the principal axis
        # Calculate the center of the object
        center = np.mean(self.position, axis=0) # axis=0 means that we take the mean of the two points in the position tuple
        # Calculate the rotation matrix to align the principal axis with the z-axis
        z_axis = principal_axis / np.linalg.norm(principal_axis)

        # Calculate the rotation matrix using the axis-angle representation
        # Calculate the angle between the z-axis and the principal axis
        sin_alpha = np.linalg.norm(np.cross(z_axis, np.array([0, 0, 1])))
        alpha = np.arcsin(sin_alpha)

        # Calculate the rotation matrix
        rotation_matrix = np.array([
            [np.cos(alpha), -np.sin(alpha), 0],
            [np.sin(alpha), np.cos(alpha), 0],
            [0, 0, 1]
        ])

        # Rotate the point to align with the new coordinate system
        rotated_point = np.dot(rotation_matrix, point - center)
        # Check if the rotated point is within the cylinder   
        r = np.linalg.norm(rotated_point[[0, 1]]) <= self.radius  # Check radial distance
        z = 0  # Check z-coordinate
        if self.position[0][2] > self.position[1][2]:
            z = self.position[0][2] <= rotated_point[2] <= self.position[1][2]
        else:
            z = self.position[1][2] <= rotated_point[2] <= self.position[0][2]
        return r and z  # Both conditions must be true to be inside the object
    

    def rotate(self, theta: float, axis: str):
        """
        Rotate the object around the center of the object and a given axis by a given angle.

        :param theta: Angle of rotation in radians.
        :param axis: Axis of rotation as a string, "x", "y", "z".
        :return: Rotated object position.
        """
        # Define the rotation matrix based on the specified axis
        if axis == 'x':
            rotation_matrix = np.array([
                [1, 0, 0],
                [0, np.cos(theta), -np.sin(theta)],
                [0, np.sin(theta), np.cos(theta)]
            ])
        elif axis == 'y':
            rotation_matrix = np.array([
                [np.cos(theta), 0, np.sin(theta)],
                [0, 1, 0],
                [-np.sin(theta), 0, np.cos(theta)]
            ])
        elif axis == 'z':
            rotation_matrix = np.array([
                [np.cos(theta), -np.sin(theta), 0],
                [np.sin(theta), np.cos(theta), 0],
                [0, 0, 1]
            ])
        else:
            raise ValueError("Axis must be 'x', 'y' or 'z'.")

        # Rotate the object position
        center = np.mean(self.position, axis=0)  # Center of the object
        # Rotate the position of the object around the center
        rotated_position = [np.dot(rotation_matrix, self.position[i] - center) + center for i in range(2)]
        # Update the position of the object
        self.position = (rotated_position[0], rotated_position[1])
        return rotated_position
    

    def draw_3D(self, ax, color='blue', alpha=0.5, label=None):
        """
        Draws the 3D representation of the object using matplotlib.

        :param ax: The matplotlib axis object to plot on.
        :param color: Color of the cylinder.
        :param alpha: Transparency of the cylinder.
        :param label: Label for the detector in the legend.
        """
        # Get the principal axis and normalize it
        principal_axis = np.array(self.position[1]) - np.array(self.position[0])
        height = np.linalg.norm(principal_axis)
        direction = principal_axis / height
        
        # Get center point of the cylinder
        center = np.mean(self.position, axis=0)
        
        # Create orthogonal vectors to the principal axis
        # First, find a vector perpendicular to the direction
        if np.abs(direction[0]) < np.abs(direction[1]) and np.abs(direction[0]) < np.abs(direction[2]):
            ortho = np.array([1, 0, 0])
        elif np.abs(direction[1]) < np.abs(direction[0]) and np.abs(direction[1]) < np.abs(direction[2]):
            ortho = np.array([0, 1, 0])
        else:
            ortho = np.array([0, 0, 1])
        
        # Create the basis vectors for the cylinder
        v1 = np.cross(direction, ortho)
        v1 = v1 / np.linalg.norm(v1)
        v2 = np.cross(direction, v1)
        
        # Create a circle in the x-y plane
        theta = np.linspace(0, 2*np.pi, 50)
        z_cylinder = np.linspace(-height/2, height/2, 20)
        Theta, Z = np.meshgrid(theta, z_cylinder)
        
        # Cylinder surface points
        X = self.radius * np.cos(Theta)
        Y = self.radius * np.sin(Theta)
        
        # Transform points to align with principal axis
        xyz = np.zeros((len(z_cylinder), len(theta), 3))
        for i in range(len(z_cylinder)):
            for j in range(len(theta)):
                # Compute the point in the local coordinate system
                p = X[i, j] * v1 + Y[i, j] * v2 + Z[i, j] * direction
                # Transform to global coordinate system
                xyz[i, j, :] = center + p
        
        # Plot the cylinder
        ax.plot_surface(xyz[:, :, 0], xyz[:, :, 1], xyz[:, :, 2], 
                        color=color, alpha=alpha, edgecolor='none')
        
        # Plot the end caps (optional)
        phi = np.linspace(0, 2*np.pi, 50)
        r = np.linspace(0, self.radius, 20)
        Phi, R = np.meshgrid(phi, r)
        
        X_cap = R * np.cos(Phi)
        Y_cap = R * np.sin(Phi)
        
        # Transform caps to align with principal axis
        for sign, endpoint in enumerate([center - direction * height/2, 
                                        center + direction * height/2]):
            xyz_cap = np.zeros((len(r), len(phi), 3))
            for i in range(len(r)):
                for j in range(len(phi)):
                    p = X_cap[i, j] * v1 + Y_cap[i, j] * v2
                    xyz_cap[i, j, :] = endpoint + p
            
            ax.plot_surface(xyz_cap[:, :, 0], xyz_cap[:, :, 1], xyz_cap[:, :, 2], 
                            color=color, alpha=alpha, edgecolor='none')       


# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# SUBCLASS DETECTOR
# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-


class Detector(Object):
    def __init__(
            self, 
            position: tuple[list[float], list[float]], 
            radius: float, # cm
            energetic_resolution: float, 
            Z: float = 49.7, # Z_eff for NaI (effective atomic number of NaI material)
            density: float = 3.67, # g/cm^3 (Density of NaI, the material of the detector)
            molar_mass: float = 149.89 # g/mol
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
        super().__init__(position=position, radius=radius, Z=Z, density=density, molar_mass=molar_mass)
        self.energetic_resolution = energetic_resolution  # Energy resolution of the detector


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
    

    def resolution(self, energy: float) -> float:
        """
        Simulates the energy resolution of the detector by adding random fluctuations.

        :param energy: The true energy of the photon (float).
        :return: The energy after applying the detector's resolution (float).
        """
        return np.random.normal(energy, energy * self.energetic_resolution)


    def detection(self, electron: p.Electron) -> float:
        """
        Simulates the energy detected by the interaction of an electron.

        :param electron: The incident electron.
        :return: The detected energy (simulated by applying resolution).
        """
        return self.resolution(electron.energy) + self.exponential_background(electron.energy)    


# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# SUBCLASS TARGET
# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-


class Target(Object): 
    def __init__(
            self, 
            position: tuple[list[float], list[float]], 
            radius: float, # cm
            Z: float = 29, # Z of Cu
            density: float = 8.96, # g/cm^3 (Density of Cu)
            molar_mass: float = 63.55 # g/mol (Molar mass of Cu)
        ):
        """
        Initialize the target with its position, dimensions, density, and atomic number.
        
        :param position: Position of the target as a 3D vector (list of floats).
        :param radius: Radius of the target in meters.
        :param width: Width (height) of the target in meters.
        :param density: Density of the target material in g/cm^3 (default: density of Pb).
        :param Z: Atomic number of the target material (default: Z of Pb).
        """
        super().__init__(position=position, radius=radius, Z=Z, density=density, molar_mass=molar_mass)

    
    def __repr__(self):
        """
        String representation for the Target object. This is called when printing the object.

        :return: A string representing the target's main properties.
        """
        return f"Target(Position={self.position}, Radius={self.radius}, Z={self.Z}, Density={self.density}, Molar mass={self.molar_mass})"
    
    
    def info(self):
        """
        Print the detailed information of the target.

        :return: None
        """ 
        print("Target Information:")
        print(f"Position: {self.position}")
        print(f"Radius: {self.radius} m")
        print(f"Z: {self.Z}")
        print(f"Density: {self.density}")
        print(f"Molar mass: {self.molar_mass}")
