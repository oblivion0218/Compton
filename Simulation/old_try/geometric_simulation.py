import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# Constants and initial parameters
#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

source_position = np.array([0, 0, 0])  # Position of the source (0,0,0)

# Detector parameters
detector_ugo_radius = 1.27  # Radius in cm (1 inch = 2.54 cm)
detector_ugo_width = 2.54  # Height in cm
detector_franco_radius = 2.54  # Radius in cm
detector_franco_width = 5.08  # Height in cm

ugo_position = np.array([0, 5, 0])  # Distance from source (15 cm along y-axis)
franco_position = np.array([0, -5, 0])  # Distance from source (-15 cm along y-axis)

#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# Photons positions
#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

# Function to generate random directions
def random_directions(n): # n = number of photons
    phi = np.random.uniform(0, 2 * np.pi, n)
    theta = np.random.uniform(0, np.pi, n)

    x = np.sin(theta) * np.cos(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(theta)

    return np.vstack((x, y, z)).T

# Function to calculate hit position
def hit_position(direction, distance, source_pos=source_position):
    x = source_pos[0] + direction[0] * distance
    y = source_pos[1] + direction[1] * distance
    z = source_pos[2] + direction[2] * distance
    return np.array([x, y, z])

# Function to check if a photon is within a detector's volume
def is_in_detector(hit_points, detector, source_pos=source_position):  
    if detector == "ugo":
        detector_pos = ugo_position
        detector_radius = detector_ugo_radius
    elif detector == "franco":
        detector_pos = franco_position
        detector_radius = detector_franco_radius

    hits = 0

    for hit_point in hit_points:     

        if detector_pos[1] * hit_point[1] > 0 and np.sqrt((hit_point[0] - source_pos[0]) ** 2 + (hit_point[2] - source_pos[2]) ** 2) <= detector_radius:
            hits += 1

    return hits

#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# Visualization
#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

# Function to draw a detector as cylinder in 2D projection
def draw_detector(ax, position, radius, width, plane='xy', color='red', label=None):

    if plane == 'xy':
        # Cylinder projection as rectangle on XY plane
        rect = patches.Rectangle(
            (position[0] - radius, position[1] if position[1] > 0 else position[1] - width),
            2 * radius,  # Larghezza del rettangolo (diametro del cilindro)
            width,      # Altezza del rettangolo (profondità del cilindro)
            color=color, alpha=0.5, label=label
        )
        ax.add_patch(rect)

    elif plane == 'xz':
        # Cylinder projection as circle on XZ plane
        circle = patches.Circle(
            (position[0], position[2]),
            radius,  # Raggio del cilindro
            color=color, alpha=0.5, label=label
        )
        ax.add_patch(circle)

    elif plane == 'yz':
        # Cylinder projection as rectangle on YZ plane
        rect = patches.Rectangle(
            (position[1] if position[1] > 0 else position[1] - width, position[2] - radius),
            width,  # Larghezza del rettangolo (diametro del cilindro)
            2 * radius,      # Altezza del rettangolo
            color=color, alpha=0.5, label=label
        )
        ax.add_patch(rect)

    else:
        raise ValueError("Plane must be 'xy', 'xz' or 'yz'.")

    # Adatta i limiti degli assi per mostrare correttamente le proiezioni
    ax.set_aspect('equal', adjustable='datalim')

def draw_cylinder(ax, base_center, radius, width, axis='z', color='r', alpha=0.5, label=None): # alpha indica la trasparenza
    theta = np.linspace(0, 2 * np.pi, 100)
    x = radius * np.cos(theta)
    y = radius * np.sin(theta)

    def linspace(base, width): 
        if base > 0:
            return np.linspace(base, base + width, 50)
        elif base < 0:
            return np.linspace(base, base - width, 50)
    
    if axis == 'z':
        z = linspace(base_center[2], width)
        X, Z = np.meshgrid(x, z)
        Y, _ = np.meshgrid(y, z)
        X += base_center[0]
        Y += base_center[1]
    elif axis == 'x':
        z = linspace(base_center[0], width)
        Z, X = np.meshgrid(x, z)
        Y, _ = np.meshgrid(y, z)
        Z += base_center[2]
        Y += base_center[1]
    elif axis == 'y':
        z =  linspace(base_center[1], width)
        Z, Y = np.meshgrid(x, z)
        X, _ = np.meshgrid(y, z)
        Z += base_center[2]
        X += base_center[0]
    else:
        raise ValueError("Axis must be 'x', 'y' or 'z'")
    
    # Disegna la superficie
    ax.plot_surface(X, Y, Z, color=color, alpha=alpha, edgecolor='none', label=label)

def norm(v, w):
    return np.sqrt((v[0] - w[0]) ** 2 + (v[1] - w[1]) ** 2 + (v[2] - w[2]) ** 2)

# Function to visualize hit positions in 3D and projections with detection count
def visualize_coincidence(fileNamePNG, directions, distance, source_pos= source_position, ugo_pos= ugo_position, ugo_radius= detector_ugo_radius, ugo_width= detector_ugo_width, franco_pos= franco_position, franco_radius= detector_franco_radius, franco_width= detector_franco_width):
    # Distances to the detectors
    ugo_distance = norm(ugo_pos, source_position)
    franco_distance = norm(franco_pos, source_position)
    
    # Calculate hit positions
    hit_points = [hit_position(dir, ugo_distance) for dir in directions]
    hit_points = np.array(hit_points)
    
    # Count photons detected by each detector
    ugo_hits = is_in_detector(hit_points, "ugo", source_pos)
    franco_hits = is_in_detector(hit_points, "franco", source_pos)
    
    # 3D Plot
    fig = plt.figure(figsize=(16, 12))
    ax = fig.add_subplot(221, projection='3d')
    ax.scatter(hit_points[:, 0], hit_points[:, 1], hit_points[:, 2], s=1, label='Hit Points')
    draw_cylinder(ax, ugo_pos, ugo_radius, ugo_width, axis='y', color='red', alpha=0.6, label='Ugo Detector')
    draw_cylinder(ax, franco_pos, franco_radius, franco_width, axis='y', color='blue', alpha=0.6, label='Franco Detector')
    ax.set_xlabel('X (cm)')
    ax.set_ylabel('Y (cm)')
    ax.set_zlabel('Z (cm)')
    ax.set_title('3D Photon Hit Positions')
    ax.legend(loc='upper right')
    
    # 2D Projections
    ax_xy = fig.add_subplot(222)
    ax_xy.scatter(hit_points[:, 0], hit_points[:, 1], s=1, label='Hits in XY')
    draw_detector(ax_xy, position=ugo_position, radius=ugo_radius, width=ugo_width, plane='xy', color='red', label='Ugo Detector')
    draw_detector(ax_xy, position=franco_position, radius=franco_radius, width=franco_width, plane='xy', color='blue', label='Franco Detector')
    ax_xy.set_xlabel('X (cm)')
    ax_xy.set_ylabel('Y (cm)')
    ax_xy.set_title('Projection on XY Plane')
    ax_xy.legend(loc='upper right')

    ax_xz = fig.add_subplot(223)
    ax_xz.scatter(hit_points[:, 0], hit_points[:, 2], s=1, label='Hits in XZ')
    draw_detector(ax_xz, position=ugo_position, radius=ugo_radius, width=ugo_width, plane='xz', color='red', label='Ugo Detector')
    draw_detector(ax_xz, position=franco_position, radius=franco_radius, width=franco_width, plane='xz', color='blue', label='Franco Detector')
    ax_xz.set_xlabel('X (cm)')
    ax_xz.set_ylabel('Z (cm)')
    ax_xz.set_title('Projection on XZ Plane')
    ax_xz.legend(loc='upper right')

    ax_yz = fig.add_subplot(224)
    ax_yz.scatter(hit_points[:, 1], hit_points[:, 2], s=1, label='Hits in YZ')
    draw_detector(ax_yz, position=ugo_position, radius=ugo_radius, width=ugo_width, plane='yz', color='red', label='Ugo Detector')
    draw_detector(ax_yz, position=franco_position, radius=franco_radius, width=franco_width, plane='yz', color='blue', label='Franco Detector')
    ax_yz.set_xlabel('Y (cm)')
    ax_yz.set_ylabel('Z (cm)')
    ax_yz.set_title('Projection on YZ Plane')
    ax_yz.legend(loc='upper right')

    # Add detection counts to the plots
    fig.suptitle(f"Photons Detected:\nUgo Detector: {ugo_hits} photons\nFranco Detector: {franco_hits} photons\nNumber of photons: {len(directions)}", fontsize=16)
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])  # Adjust for title
    plt.savefig(fileNamePNG)
