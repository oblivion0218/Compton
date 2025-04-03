import numpy as np
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
from lib import detector as d
from lib import particles as p
from lib import experiments as e
from lib import source as s
from lib import visualization as v

# File path to save the output spectrum plot
file_path = "/mnt/c/Users/User/Desktop/info/Compton/Simulation/plots/"

# Initialize detectors with their respective positions and dimensions
ugo = d.Detector(([0, -16, 0], [0, -18.58, 0]), 1.27, 0.0903)  # Detector "Ugo"
franco = d.Detector(([0, 30.5, 0], [0, 35.68, 0]), 2.54, 0.0695)  # Detector "Franco"

number_of_photons = 1000000  # Number of photons to simulate

# Simulate spectroscopy measurements and plot the energy spectrum 
energies = e.spectroscopy_measurement(number_of_photons, ugo, s.Source(), True)
v.plot_energy_spectrum(energies, file_path + "spectrum/spettroscopy_ugo.png", 10000)
energies = e.spectroscopy_measurement(number_of_photons, franco, s.Source(), True)
v.plot_energy_spectrum(energies, file_path + "spectrum/spettroscopy_franco.png", 10000)
print("Spectroscopy measurements completed.")

# Simulate photon emission from a source
source = s.Source()  # Create a source object
photons = source.photon_emission(1000, np.arctan(1.27/16), 2 * np.pi)  # Simulate photon emission
[photon.propagation(16) for photon in photons]  # Propagate the photons

franco.rotate((5/18) * np.pi, [0, 5.5, 0], "z")

# Visualize the 3D photon hit positions and detectors for coincidence measurement
v.visualization_3D(file_path + "3D_plots/3D_visualization.png", [ugo, franco], photons, source)
print("3D visualization completed.")

target = d.Target(([0, 5, 0], [0, 6, 0]), 3)  # Create a target object

# # Generate interactive 3D visualization with Plotly
interactive_html_path = file_path + "3D_plots/interactive_3D_visualization.html"
v.visualization_3D_plotly(interactive_html_path, [ugo, franco], photons, source, target)
print("Interactive 3D visualization completed.")

print("End")

