import numpy as np
from . import particles as p
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


    def principal_axis(self) -> np.ndarray:
        """
        Calculate the principal axis of the object, which is the vector connecting the two points in the position tuple.

        :return: The principal axis as a numpy array.
        """
        return np.array(self.position[1]) - np.array(self.position[0])
    

    def center(self) -> np.ndarray:
        """
        Calculate the center of the object.

        :return: The center of the object as a numpy array.
        """
        return np.mean(self.position, axis=0) # axis=0 means that we take the mean of the two points in the position tuple
    

    # def is_inside(self, point: np.ndarray) -> bool:
    #     """
    #     Checks if a given point (position of a photon) is inside the object.

    #     :param point: The position of the particle (numpy array).
    #     :return: True if the point is within the object, False otherwise.
    #     """
    #     # Ensure point is a numpy array
    #     point = np.array(point)
        
    #     # Define the principal axis of the object
    #     principal_axis = self.principal_axis()
    #     axis_length = np.linalg.norm(principal_axis)
        
    #     # Get center and normalize axis
    #     center = self.center()
    #     if axis_length == 0:
    #         return False  # Degenerate cylinder
        
    #     axis_unit = principal_axis / axis_length
        
    #     # Calculate the rotation matrix to align the principal axis with the y-axis
    #     # Note: The code was using y_axis but trying to align with different axes
        
    #     # Calculate angle between principal axis and y-axis
    #     cos_alpha = np.linalg.norm(np.dot(axis_unit, np.array([0, 1, 0])))
    #     alpha = np.arcsin(cos_alpha)
        
    #     # Calculate rotation matrix
    #     rotation_matrix = np.array([
    #         [np.cos(alpha), np.sin(alpha), 0],
    #         [-np.sin(alpha), np.cos(alpha), 0], 
    #         [0, 0, 1]
    #     ])
        
    #     # Rotate the point to align with the rotated coordinate system
    #     rotated_point = np.dot(rotation_matrix, point - center)
        
    #     # Check radial distance (in x-z plane of rotated coordinates)
    #     r = np.linalg.norm(rotated_point[[0, 2]]) <= self.radius
        
    #     # Check if point is within cylinder height
    #     # Half-length of cylinder along principal axis
    #     half_length = axis_length / 2
        
    #     # Check if y-coordinate in rotated system is within half-length
    #     h = abs(rotated_point[1]) <= half_length
        
    #     return r and h  # Both conditions must be true


    def is_inside(self, point: np.ndarray) -> bool:
        """
        Checks if a given point is inside the cylindrical object defined by two endpoints and a radius.

        :param point: The position to check (numpy array).
        :return: True if the point is inside the cylinder, False otherwise.
        """
        axis = self.principal_axis()
        axis_length = np.linalg.norm(axis)
        if axis_length == 0:
            return False  # degenerate cylinder
        axis_unit = axis / axis_length # Normalized direction of the cylinder axis
         
        v = point - self.position[0]  # Vector from one end of the cylinder to the point
        
        h = np.dot(v, axis_unit) # Projection of v onto the axis
        # Check if the projection is within the cylinder length
        if h < 0 or h > axis_length:
            return False  

        closest_point_on_axis = self.position[0] + h * axis_unit # Closest point on the axis to the point

        radial_vector = point - closest_point_on_axis # Vector from the closest point on the axis to the point
        radial_distance = np.linalg.norm(radial_vector) # Radial distance from the axis to the point

        return radial_distance <= self.radius
    

    def rotate(self, theta: float, rotation_center: list[float], axis: str):
        """
        Rotate the object around the center of the object and a given axis by a given angle.

        :param theta: Angle of rotation in radians.
        :param rotation_center: Center of rotation as a list of floats.
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

        rotation_center = np.array(rotation_center)

        # Rotate the position of the object around the center
        rotated_position = [np.dot(rotation_matrix, self.position[i] - rotation_center) + rotation_center for i in range(2)]
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
        principal_axis = self.principal_axis()
        height = np.linalg.norm(principal_axis)
        direction = principal_axis / height
        
        # Get center point of the cylinder
        center = self.center()
        
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


    def draw_plotly_3D(self, color='blue', alpha=0.7, name=None):
        """
        Creates a Plotly mesh3d representation of the cylindrical object.
        
        Args:
            color: Color of the cylinder
            alpha: Transparency of the cylinder (0-1)
            name: Name for the trace in the legend
            
        Returns:
            Plotly mesh3d trace object
        """
        import plotly.graph_objects as go
        
        # Get the principal axis and normalize it
        principal_axis = self.principal_axis()
        height = np.linalg.norm(principal_axis)
        if height == 0:
            return None  # Can't create a cylinder with zero height
            
        direction = principal_axis / height
        
        # Get center point of the cylinder
        center = self.center()
        
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
        
        # Number of points around the cylinder circumference
        n_theta = 36
        
        # Create points for the vertices of the mesh
        theta = np.linspace(0, 2*np.pi, n_theta)
        vertices = []
        # Create vertices for top and bottom circles
        for t in theta:
            circle_pt = self.radius * (v1 * np.cos(t) + v2 * np.sin(t))
            # Bottom circle point
            vertices.append(list(center - direction * height/2 + circle_pt))
            # Top circle point
            vertices.append(list(center + direction * height/2 + circle_pt))
        
        # Create indices for triangular faces
        i_vals = []
        j_vals = []
        k_vals = []
        
        # Connect vertices to form triangles
        for i in range(n_theta):
            i0 = 2 * i
            i1 = 2 * i + 1
            i2 = 2 * ((i + 1) % n_theta)
            i3 = 2 * ((i + 1) % n_theta) + 1
            
            # Triangle 1 (connecting points on adjacent sides)
            i_vals.append(i0)
            j_vals.append(i2)
            k_vals.append(i1)
            
            # Triangle 2 (for the other half of the quad)
            i_vals.append(i1)
            j_vals.append(i2)
            k_vals.append(i3)
        
        # Add triangles for the bottom cap
        # First add center of bottom cap as a vertex
        bottom_center_idx = len(vertices)
        vertices.append(list(center - direction * height/2))
        
        for i in range(n_theta):
            i0 = 2 * i
            i2 = 2 * ((i + 1) % n_theta)
            
            i_vals.append(bottom_center_idx)
            j_vals.append(i0)
            k_vals.append(i2)
        
        # Add triangles for the top cap
        # Add center of top cap as a vertex
        top_center_idx = len(vertices)
        vertices.append(list(center + direction * height/2))
        
        for i in range(n_theta):
            i1 = 2 * i + 1
            i3 = 2 * ((i + 1) % n_theta) + 1
            
            i_vals.append(top_center_idx)
            j_vals.append(i3)
            k_vals.append(i1)
        
        # Create the mesh3d trace
        return go.Mesh3d(
            x=[v[0] for v in vertices],
            y=[v[1] for v in vertices],
            z=[v[2] for v in vertices],
            i=i_vals,
            j=j_vals,
            k=k_vals,
            color=color,
            opacity=alpha,
            name=name or "Cylinder"
        )


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
        return np.random.normal(energy, energy * (self.energetic_resolution / 2.355))


    def detection(self, electron: p.Electron) -> float:
        """
        Simulates the energy detected by the interaction of an electron.

        :param electron: The incident electron.
        :return: The detected energy (simulated by applying resolution).
        """
        return self.resolution(electron.energy)


# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# SUBCLASS TARGET
# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-


class Target(Object): 
    def __init__(
            self, 
            position: tuple[list[float], list[float]], 
            radius: float, # cm
            Z: float = 29, # Z of Cu
            density: float = 8.935, # g/cm^3 (Density of Cu)
            molar_mass: float = 63.546 # g/mol (Molar mass of Cu)
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
