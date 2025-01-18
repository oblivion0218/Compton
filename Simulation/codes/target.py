import numpy as np
import random
import particles as p
import matplotlib.patches as patches

class Target:
    def __init__(
            self, 
            position: list[float], 
            radius: float, # cm
            width: float, #cm
            Z: float, 
            density: float, # g/cm^3
            molar_mass: float # g/mol
        ):
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
        self.Z = Z  # Atomic number of the detector 
        self.density = density
        self.molar_mass = molar_mass

    def __repr__(self):
        """
        String representation for the Target object.
        """
        return f"Target(Position={self.position}, Radius={self.radius}, Width={self.width}, Z={self.Z}, Density={self.density}, Molar mass={self.molar_mass})"
    
    def info(self):
        """
        Print the detailed information of the target.
        """
        print("Target Information:")
        print(f"Position: {self.position}")
        print(f"Radius: {self.radius} m")
        print(f"Width: {self.width} m")
        print(f"Z: {self.Z}")
        print(f"Density: {self.density}")
        print(f"Molar mass: {self.molar_mass}")
    
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