import numpy as np
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
from lib import detector as d
from lib import particles as p
from lib import experiments as e
from lib import source as s

# File path to save the output spectrum plot
file_path = "/mnt/c/Users/User/Desktop/info/Compton/Simulation/plots/"


def plot_energy_spectrum(energies: list[float], fileNamePNG: str, bins: int=100, title: str="Energy Spectrum"):
    """
    Plots a histogram of the detected photon energies.
    
    :param energies: List of detected photon energies in keV.
    :param fileNamePNG: File path to save the plot as a PNG image.
    :param bins: Number of bins for the histogram.
    :param title: Title of the plot.
    """
    # Filter out energies less than or equal to 1 keV, which might be noise or irrelevant
    energies_def = [energy for energy in energies if energy > 1]

    if not energies_def:
        print("No energies to plot.")
        return

    # Create a histogram of detected energies
    plt.figure(figsize=(10, 6))
    plt.hist(energies_def, bins=bins, color="blue", alpha=0.7, edgecolor="blue")
    
    # PP ---> Photopeak, CE ---> Compton Edge
    # Mark the photopeak (PP) and Compton Edge (CE) for the energies of interest
    plt.axvline(511, color='red', linestyle='--', linewidth=1.5, label='PP for 511 keV')
    plt.axvline(1274, color='red', linestyle='--', linewidth=1.5, label='PP for 1274 keV')
    plt.axvline(341, color='orange', linestyle='--', linewidth=1.5, label='CE for 511 keV') # 2/3 511
    plt.axvline(1061, color='orange', linestyle='--', linewidth=1.5, label='CE for 1274 keV') 
    plt.legend(loc="upper right")

    # Set plot title and labels
    plt.title(title, fontsize=16)
    plt.xlabel("Energy (keV)", fontsize=14)
    plt.ylabel("Counts", fontsize=14)
    plt.yscale("log")  # Log scale for better visibility of data
    plt.grid(alpha=0.4)  # Grid for readability
    plt.tight_layout()  # Adjust layout to fit everything
    plt.savefig(fileNamePNG)  # Save the plot as PNG
    plt.close()  # Close the plot to free up resources


def visualization_3D(fileNamePNG: str, detectors: d.Detector, photons: list[p.Photon], source: s.Source, target: d.Target=None):
    """
    Creates a 3D visualization of the photon hit points and detector positions.
    
    :param fileNamePNG: Path to save the output 3D visualization image.
    :param detectors: List of detector objects.
    :param photons: List of photon objects, each with a position attribute.
    :param target: Optional target object to visualize.
    """
    hit_points = np.array([photon.position for photon in photons])  # Extract photon hit positions
    
    # Create figure with a single 3D plot
    fig = plt.figure(figsize=(14, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    # Plot hit points
    ax.scatter(hit_points[:, 0], hit_points[:, 1], hit_points[:, 2], 
               s=10, color='red', marker='x', label='Hit Points')
    
    # Draw detectors and target
    [detector.draw_3D(ax) for detector in detectors]
    if target:
        target.draw_3D(ax)

    # Draw source
    source.draw_3D(ax)
    
    # Add labels and title
    ax.set_xlabel('X (cm)', fontsize=12)
    ax.set_ylabel('Y (cm)', fontsize=12)
    ax.set_zlabel('Z (cm)', fontsize=12)
    ax.set_title('3D Photon Hit Positions', fontsize=14)
    
    # Improve visibility
    ax.legend(loc='upper right', fontsize=10)
    ax.grid(True, alpha=0.3)
    
    # Set equal aspect ratio for better visualization
    ax.set_box_aspect([1, 1, 1])
    
    # Adjust layout and save
    plt.tight_layout()
    plt.savefig(fileNamePNG, dpi=300)
    plt.close()


def visualization_3D_plotly(fileNameHTML: str, detectors: d.Detector, photons: list[p.Photon], source: s.Source=None, target: d.Target=None):
    """
    Creates an interactive 3D visualization of the photon hit points and detector positions using Plotly.
    
    :param fileNameHTML: Path to save the output 3D visualization as HTML.
    :param detectors: List of detector objects.
    :param photons: List of photon objects, each with a position attribute.
    :param source: Optional source object to visualize.
    :param target: Optional target object to visualize.
    :return: The Plotly figure object (can be displayed in notebooks)
    """
    
    # Extract position data from photons
    photon_x = [photon.position[0] for photon in photons]
    photon_y = [photon.position[1] for photon in photons]
    photon_z = [photon.position[2] for photon in photons]
    
    # Create the figure
    fig = go.Figure()
    
    # Add photon hit points
    fig.add_trace(go.Scatter3d(
        x=photon_x, y=photon_y, z=photon_z,
        mode='markers',
        marker=dict(size=3, color='red', symbol='cross'),
        name='Photon Hit Points'
    ))
    
    # Add detectors using the draw_plotly_3D method
    for i, detector in enumerate(detectors):
        # Alternate colors for different detectors
        color = 'blue' if i % 2 == 0 else 'green'
        name = f"Detector {i+1}"
        
        # Get mesh trace from detector and add to figure
        cylinder_trace = detector.draw_plotly_3D(color=color, name=name)
        if cylinder_trace is not None:  # Check if the trace is valid
            fig.add_trace(cylinder_trace)
    
    # Add source if provided
    if source:
        # Get source trace and add to figure
        source_trace = source.draw_plotly_3D(color='black', name='Source')
        if source_trace is not None:  # Check if the trace is valid
            fig.add_trace(source_trace)
        
    # Add target if provided
    if target:
        target_trace = target.draw_plotly_3D(color='gray', name='Target')
        if target_trace is not None:  # Check if the trace is valid
            fig.add_trace(target_trace)
    
    # Update layout
    fig.update_layout(
        title='3D Photon Hit Positions',
        scene=dict(
            xaxis_title='X (cm)',
            yaxis_title='Y (cm)',
            zaxis_title='Z (cm)',
            aspectmode='data'  # Ensure proper aspect ratio
        ),
        legend=dict(
            yanchor="top",
            y=0.99,
            xanchor="right",
            x=0.99
        )
    )
    
    # Save as HTML for interactive viewing
    pio.write_html(fig, file=fileNameHTML)
    print(f"Interactive 3D visualization saved to {fileNameHTML}")
    
    return fig
