import pandas as pd
import matplotlib.pyplot as plt
import numpy as np # Import numpy for CDF calculation
import os
from tqdm import tqdm

# Define the base path and create plot directory if it doesn't exist
file_path = "/mnt/c/Users/User/Desktop/info/Compton/Simulations/Python_simulation/interaction_point/"


angles = [35, 40, 50, 60, 70, 80, 90, 100, 110, 120]
N_photons = 1000000 # Assuming this is a constant for normalization across simulations
data = {35: [], 40: [], 50: [], 60: [], 70: [], 80: [], 90: [], 100: [], 110: [], 120: []}
Omega_orizontal = {35: [], 40: [], 50: [], 60: [], 70: [], 80: [], 90: [], 100: [], 110: [], 120: []}
Omega_vertical = {35: [], 40: [], 50: [], 60: [], 70: [], 80: [], 90: [], 100: [], 110: [], 120: []}
Omega_mean = {35: [], 40: [], 50: [], 60: [], 70: [], 80: [], 90: [], 100: [], 110: [], 120: []}

for filename in os.listdir(file_path + "data/"):
    # Extract angle string: e.g., "lambda_simulation_30.txt" -> "30"
    angle_str = filename[len("IP_simulation_"):-len(".txt")]
    angle = int(angle_str)
    
    with open(file_path + "data/" + filename, 'r') as f:
        ip = []
        for i, line in enumerate(f):
            # Clean the line by removing brackets and splitting by commas or spaces
            clean_line = line.strip().replace('[', '').replace(']', '')
            
            # Try different splitting methods since we don't know the exact format
            if ',' in clean_line:
                # If commas are present, split by commas
                values = [float(val.strip()) for val in clean_line.split(',')]
            else:
                # Otherwise split by whitespace
                values = [float(val) for val in clean_line.split()]
            
            # Make sure we have 3 values (x, y, z)
            if len(values) == 3:
                ip.append(values)
            else:
                print(f"Warning: Line {i+1} in file {filename} does not contain exactly 3 values: {line}")
        
        data[angle] = ip

target_center = [0, 5.5, 0]
detector_center_0 = [0, 33.04, 0]
d_target_detector = np.abs(target_center[1] - detector_center_0[1])
detector_radius = 2.54

with open("solid_angle.txt", "a") as f:
    f.write("Angle\tMean\tMedian\tIQR\n")

for angle in tqdm(angles, desc="Calculating Omega", unit="angle"):
    detector_center_rotated = [target_center[0] - d_target_detector * np.sin(angle), target_center[1] + d_target_detector * np.cos(angle), 0]

    detector_extreme_orizontal = np.array([[detector_center_rotated[0] + detector_radius * np.cos(angle), detector_center_rotated[1] + detector_radius * np.sin(angle), detector_center_rotated[2]],
                                    [detector_center_rotated[0] - detector_radius * np.cos(angle), detector_center_rotated[1] - detector_radius * np.sin(angle), detector_center_rotated[2]]])

    detector_extreme_vertical = np.array([[detector_center_rotated[0], detector_center_rotated[1], detector_center_rotated[2] + detector_radius],
                                    [detector_center_rotated[0], detector_center_rotated[1], detector_center_rotated[2] - detector_radius]])

    for d in data[angle]:
        ip = np.array([d[0], d[1], d[2]])

        A_orizontal = np.linalg.norm(detector_extreme_orizontal[0] - ip)
        B_orizontal = np.linalg.norm(detector_extreme_orizontal[1] - ip)

        A_vertical = np.linalg.norm(detector_extreme_vertical[0] - ip)
        B_vertical = np.linalg.norm(detector_extreme_vertical[1] - ip)

        C = 2 * detector_radius

        beta_orizontal = np.arccos((A_orizontal ** 2 + B_orizontal ** 2 - C ** 2) / (2 * A_orizontal * B_orizontal))
        beta_vertical = np.arccos((A_vertical ** 2 + B_vertical ** 2 - C ** 2) / (2 * A_vertical * B_vertical))

        omega_orizontal = 2 * np.pi * (1 - np.cos(beta_orizontal / 2))
        omega_vertical = 2 * np.pi * (1 - np.cos(beta_vertical / 2))

        Omega_orizontal[angle].append(omega_orizontal)
        Omega_vertical[angle].append(omega_vertical)
        Omega_mean[angle].append((omega_orizontal + omega_vertical) / 2)


    # print an histogram of the Omega_orizontal and Omega_vertical for each angle
    plt.figure(figsize=(12, 6))
    plt.hist(Omega_orizontal[angle], bins=100, alpha=0.5, color='blue', label='Omega_orizontal', density=True)
    plt.hist(Omega_vertical[angle], bins=100, alpha=0.5, color='red', label='Omega_vertical', density=True)
    plt.hist(Omega_mean[angle], bins=100, alpha=0.5, color='green', label='Omega_mean', density=True)

    mean_omega = np.mean(Omega_mean[angle])
    median_omega = np.median(Omega_mean[angle])

    # Calculate IQR instead of STD
    q75 = np.percentile(Omega_mean[angle], 75)
    q25 = np.percentile(Omega_mean[angle], 25)
    iqr_omega = q75 - q25

    plt.text(0.1, 0.9, f'Mean: {mean_omega:.5f}\nMedian: {median_omega:.5f}\nIQR: {iqr_omega:.5f}', transform=plt.gca().transAxes, fontsize=12,
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))

    print(angle, mean_omega, median_omega, iqr_omega)

    plt.xlabel('Omega')
    plt.ylabel('Probability Density')  # Changed y-label to reflect normalization
    plt.title(f'Normalized Omega Histogram for Angle {angle}')  # Updated title
    plt.legend()
    plt.grid()
    plt.savefig(file_path + "plots/Omega_histogram_" + str(angle) + ".png")
    plt.close()

