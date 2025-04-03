import numpy as np
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go
from . import particles as p

# Define energies and source properties
energies = {511: 0.903, 1274: 0.097}  # Energies in keV with corresponding probabilities
source_position = np.array([0, 0, 0])  # Source position in 3D space
activity = 127000  # Activity of the source in Bq (decays per second)

class Source:        
    def __init__(self, energies: dict = energies, position: list[float] = source_position, activity: int = activity):
        """
        Initialize the source with energy levels, position, and activity.
        
        :param energies: Dictionary of energy levels with their probabilities.
                         Example: {511: 0.903, 1274: 0.097}
        :param position: Position of the source as a 3D vector.
        :param activity: Activity of the source in Bq.
        """
        self.energies = energies
        self.position = np.array(position)
        self.activity = activity


    def info(self):
        """
        Print the information about the source: energies, position, and activity.
        """
        print("Possible energies:")
        energies = self.energies.keys()
        for energy in energies:
            print(f"{energy} keV with probability of {self.energies[energy] * 100} %")
        print(f"Position: {self.position}")
        print(f"Activity: {self.activity}")


    def random_energies(self, number_of_photons: int = 1) -> np.ndarray:
        """
        Generate a list of random photon energies based on the energy probabilities of the source.
        
        :param number_of_photons: The number of photons to generate.
        :return: Array of randomly chosen photon energies.
        """
        possible_energies = list(self.energies.keys())
        probabilities = [self.energies[energy] for energy in possible_energies]
        return np.random.choice(possible_energies, size=number_of_photons, p=probabilities)


    def random_directions(self, number_of_photons: int = 1) -> np.ndarray:
        """
        Generate a list of random photon directions (unit vectors).
        
        :param number_of_photons: The number of photons to generate.
        :return: Array of photon directions as unit vectors.
        """
        phi = np.random.uniform(0, 2 * np.pi, number_of_photons)  # Random azimuthal angle
        theta = np.random.uniform(0, np.pi, number_of_photons)  # Random polar angle

        # Calculate direction components (unit vectors)
        directions = np.vstack([
            np.sin(theta) * np.cos(phi),
            np.sin(theta) * np.sin(phi),
            np.cos(theta)
        ]).T
        return directions
    
    def cone_random_directions(self, phi_max: float, theta_max: float, axis: str, number_of_photons: int = 1, forward_backward: bool = True) -> np.ndarray:
        """
        Generate a list of photon directions given a range for theta and phi (unit vectors).
        
        :param phi_max: Maximum azimuthal angle in radians.
        :param theta_max: Maximum polar angle in radians.
        :param axis: Axis of photon generation ('x', 'y', or 'z').
        :param number_of_photons: The number of photons to generate.
        :param forward_backward: If True, generate photons in both forward and backward directions.
        """
        phi = np.random.uniform(0, phi_max, number_of_photons)  # Random azimuthal angle
        theta = np.random.uniform(0, theta_max, number_of_photons)  # Random polar angle

        # Calculate direction components (unit vectors)
        if axis == 'x':
            directions = np.vstack([
                np.cos(theta),
                np.sin(theta) * np.cos(phi),
                np.sin(theta) * np.sin(phi)
            ]).T
        elif axis == 'y':
            directions = np.vstack([
                np.sin(theta) * np.cos(phi),
                np.cos(theta),
                np.sin(theta) * np.sin(phi)
            ]).T
        elif axis == 'z':
            directions = np.vstack([
                np.sin(theta) * np.cos(phi),
                np.sin(theta) * np.sin(phi),
                np.cos(theta)
            ]).T
            
        if forward_backward:
            # Generate one random choice per direction
            F_or_B = np.random.choice([-1, 1], number_of_photons)
            # Apply the random choice to each direction vector
            directions = directions * F_or_B[:, np.newaxis]
           
        return directions


    def photon_emission(self, number_of_photons: int = 1, theta_max: float = None, phi_max: float = None, axis: str="y", forward_backward: bool = True) -> list:
        """
        Generate a list of Photon objects with random energies and directions.
        
        :param number_of_photons: The number of photons to generate.
        :param theta_max: Maximum polar angle in radians.
        :param phi_max: Maximum azimuthal angle in radians.
        :param axis: Axis of photon generation ('x', 'y', or 'z').
        :param forward_backward: If True, generate photons in both forward and backward directions.
        """
        energies = self.random_energies(number_of_photons)  # Generate random photon energies
        if theta_max is not None and phi_max is not None:
            directions = self.cone_random_directions(phi_max, theta_max, axis, number_of_photons, forward_backward)
        else:
            directions = self.random_directions(number_of_photons)
            
        return [p.Photon(energy, direction) for energy, direction in zip(energies, directions)]  # Create and return Photon objects

    
    def testing_photons(self, number_of_photons: int = 1, direction: list = [0, 1, 0]) -> list:
        """
        Generate a list of Photon objects with fixed directions and random energies.
        
        :param number_of_photons: The number of photons to generate.
        :param direction: Fixed direction for all photons.
        :return: List of Photon objects.
        """
        directions_po = [np.array(direction) for _ in range(number_of_photons)]  # Create list of fixed directions
        directions = np.array(directions_po)  # Convert to NumPy array for efficient processing
        
        energies = self.random_energies(number_of_photons)  # Generate random photon energies
        return [p.Photon(energy, direction) for energy, direction in zip(energies, directions)]  # Return list of Photon objects


    def draw_3D(self, ax):
        """
        Draw the source position in 3D space.
        
        :param ax: Matplotlib 3D axis object.
        """
        ax.scatter(self.position[0], self.position[1], self.position[2], 
                   color='black', marker='o', label='Source Position')
        ax.text(self.position[0], self.position[1], self.position[2], 
                'Source', color='black', fontsize=12)


    def draw_plotly_3D(self, color='black', alpha=1.0, name=None):
        """
        Creates a Plotly representation of the source.
        
        Args:
            color: Color of the source point
            alpha: Opacity of the source marker (0-1)
            name: Name for the trace
            
        Returns:
            Plotly scatter3d trace object
        """        
        return go.Scatter3d(
            x=[self.position[0]],
            y=[self.position[1]],
            z=[self.position[2]],
            mode='markers',
            marker=dict(size=10, color=color, symbol='circle', opacity=alpha),
            name=name or "Source"
        )