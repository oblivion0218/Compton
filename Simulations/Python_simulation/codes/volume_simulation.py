import numpy as np
import matplotlib.pyplot as plt
import random
from tqdm import tqdm
from lib import interactions as i
from lib import detector as d
from lib import visualization as v


# Object initialization
# *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
angle = 40

angle_rad = angle * np.pi / 180  
target_angle_rad = - ((180 - angle)/2) * np.pi / 180  

target = d.Target(([0, -0.5, 0], [0, 0.5, 0]), 3)
target.rotate(target_angle_rad, [0, 0, 0], "z")  
v.visualization_3D_plotly("Volume_simulation.html", [], [], target=target)  # Create an empty plot for the target