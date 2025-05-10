import pandas as pd
import matplotlib.pyplot as plt
import numpy as np # Import numpy for CDF calculation
import os

# Define the base path and create plot directory if it doesn't exist
file_path = "/mnt/c/Users/User/Desktop/info/Compton/Simulations/Python_simulation/mean_free_path/"


angles = [30, 40, 50, 60, 70, 80, 90, 100, 110, 120]
N_photons = 1000000 # Assuming this is a constant for normalization across simulations
data = {30: [], 40: [], 50: [], 60: [], 70: [], 80: [], 90: [], 100: [], 110: [], 120: []}

for filename in os.listdir(file_path + "data/"):
    # Extract angle string: e.g., "lambda_simulation_30.txt" -> "30"
    angle_str = filename[len("lambda_simulation_"):-len(".txt")]
    angle = int(angle_str)
    
    with open(file_path + "data/" + filename, 'r') as f:
        lam = []
        for i, line in enumerate(f):
            lambda_val = float(line.strip())
            lam.append(lambda_val)
        data[angle] = lam

# Plot a histogram and CDF for each angle in the predefined list
for angle in angles:
    lam = data[angle]

    counts, bin_edges = np.histogram(lam, bins=100, weights=np.ones_like(lam) / N_photons)
    bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    first_bin_height = counts[0]
    line_height = first_bin_height * np.exp(-1)  # 1/3 of the first bin height

    # Plot
    plt.figure(figsize=(10, 5))
    weights = np.ones_like(lam) / N_photons
    plt.hist(lam, bins=100, alpha=0.5, color='blue', weights=weights)
    plt.plot(bin_centers, counts, color='blue', label="Mean Free Path")

    plt.axvline(1 / np.cos((np.pi - np.radians(angle)) / 2), color='red', linestyle='dashed', linewidth=1, label=f"Expected D' for D=1cm:{1 / np.cos((np.pi - np.radians(angle)) / 2):.4f} cm")

    plt.axhline(line_height, color='green', linestyle='dotted', linewidth=1, label="1/3 First Bin Height")

    # print the x value of the interseption between the line and the histogram
    x_intersection = bin_centers[np.where(np.abs(counts - line_height) == np.min(np.abs(counts - line_height)))[0][0]]
    print(x_intersection)
    plt.axvline(x_intersection, color='black', linewidth=1, label=f"Simulated lambda:{x_intersection:.4f} cm")

    plt.title(f"Histogram of Mean Free Path for {angle} degrees")
    plt.xlabel("Mean Free Path (cm)")
    plt.ylabel("Density")
    plt.legend()
    plt.savefig(file_path + "plots/histogram_" + str(angle) + ".png")
    plt.close()

