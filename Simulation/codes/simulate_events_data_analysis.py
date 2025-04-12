import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from tqdm import tqdm
from iminuit import Minuit
from iminuit.cost import LeastSquares

N_tot = 1000000  # Total number of photons in each simulation

# file_path = "/mnt/c/Users/User/Desktop/info/Compton/Simulation/simulated_events/"
file_path = "/mnt/c/Users/User/Desktop/info/Gamma-simulation/simulated_events/40_deg/"
angle = 40  # Angle in degrees
angle_rad = angle * np.pi / 180  # Convert to radians

sim_runs = []
photons_left_target = []
interacting_photons = []
photoelectric_interactions = []
compton_1st = []
compton_2nd = []
compton_3rd = []
compton_4th = []
photons_reached_detector = []
photons_observed = []
energies = []
true_energies = []

file_names = []

if os.path.isdir(file_path + "data/"):
    # Iterate through the items in the directory
    for file_name in os.listdir(file_path + "data/"):
        # Create the full file path
        full_file_path = os.path.join(file_path + "data/", file_name)
        # Add the file name if it is a file (not a directory)
        if os.path.isfile(full_file_path):
            file_names.append(file_path + "data/" + file_name)

# Process each file
for file_name in file_names:
    # Extract run number from filename
    run_number = int(os.path.basename(file_name).split("_")[0])
    sim_runs.append(run_number)
    
    # Parse the metrics from the file
    with open(file_name, 'r') as f:
        lines = f.readlines()
        
        # Extract metrics
        photons_left_target.append(int(lines[0].split(": ")[1]))
        interacting_photons.append(int(lines[1].split(": ")[1]))
        photoelectric_interactions.append(int(lines[5].split(": ")[1]))
        compton_1st.append(int(lines[9].split(": ")[1]))
        compton_2nd.append(int(lines[10].split(": ")[1]))
        compton_3rd.append(int(lines[11].split(": ")[1]))
        compton_4th.append(int(lines[12].split(": ")[1]))
        photons_reached_detector.append(int(lines[13].split(": ")[1]))
        photons_observed.append(int(lines[14].split(": ")[1]))
        
        capture_true_energies = False
        for i in range(16, len(lines)):
            line = lines[i].strip()
            if line == "True Detected Energies:":
                capture_true_energies = True
                continue
            if line and not capture_true_energies:  # Detected energies
                energies.append(float(line))
            elif line and capture_true_energies:  # True energies
                true_energies.append(float(line))

# Convert lists to numpy arrays for easier manipulation
sim_runs = np.array(sim_runs)
photons_left_target = np.array(photons_left_target)
interacting_photons = np.array(interacting_photons)
photoelectric_interactions = np.array(photoelectric_interactions)
compton_1st = np.array(compton_1st)
compton_2nd = np.array(compton_2nd)
compton_3rd = np.array(compton_3rd)
compton_4th = np.array(compton_4th)
photons_reached_detector = np.array(photons_reached_detector)
photons_observed = np.array(photons_observed)
energies = np.array(energies)
true_energies = np.array(true_energies)

plt.figure(figsize=(12, 8))

# Plotting all metrics on the same plot
plt.plot(sim_runs, interacting_photons/N_tot, 'o', color='red', label='Photon Interactions')
plt.plot(sim_runs, photoelectric_interactions/N_tot, 'o', color='green', label='Photoelectric Interactions')
plt.plot(sim_runs, compton_1st/N_tot, 'o', color='purple', label='Compton 1st')
plt.plot(sim_runs, compton_2nd/N_tot, 'o', color='orange', label='Compton 2nd')
plt.plot(sim_runs, compton_3rd/N_tot, 'o', color='brown', label='Compton 3rd')
plt.plot(sim_runs, compton_4th/N_tot, 'o', color='pink', label='Compton 4th')

plt.title('Photon Interaction probability across simulation runs')
plt.xlabel('Simulation Run')
plt.ylabel('Probability / Scaled Efficiency')
plt.legend(fontsize=12)
plt.grid(True)
plt.savefig(file_path + "plots/interaction_probabilities.png")

plt.figure(figsize=(12, 8))
plt.plot(sim_runs, photons_reached_detector/N_tot, 'o', color='red', label='Photons Reaching Detector')
plt.plot(sim_runs, photons_observed/N_tot, 'o', color='blue', label='Photons Observed')
plt.title('Photon interaction across first simulation runs')
plt.xlabel('Simulation run')
plt.ylabel('Probability / Scaled Efficiency')
plt.legend(fontsize=12)
plt.grid(True)
plt.savefig(file_path + "plots/detector_parameter.png")

detection_efficiency = photons_observed / photons_reached_detector
sorted_indices = np.argsort(sim_runs)
sim_runs_sorted = sim_runs[sorted_indices]
detection_efficiency_sorted = detection_efficiency[sorted_indices]

plt.figure(figsize=(12, 8))
plt.plot(sim_runs, detection_efficiency, 'o', color='black', label='Detection Efficiency (scaled)')
plt.title('Detection Efficiency across simulation runs')
plt.xlabel('Simulation Run')
plt.ylabel('Detection Efficiency (scaled)')
plt.grid(True)
plt.savefig(file_path + "plots/detection_efficiency.png")

# Create a new combined histogram with Freedman-Diaconis binning
plt.figure(figsize=(12, 8))

# Create histograms
plt.hist(energies, bins=100, alpha=0.4, color='blue', 
         edgecolor='black', label='Detected Energies')
plt.hist(true_energies, bins=100, color='red', 
              histtype="step", label='True Energies')

plt.title('Probability Density of Photon Energies')
plt.xlabel('Energy (keV)')
plt.ylabel('Counts')
plt.legend(fontsize=12)
plt.grid(True)
plt.savefig(file_path + "plots/photon_energies.png")

# Fit the histogram of detected energies
def gaussian(x, a, x0, sigma):
    return a * np.exp(-0.5 * ((x - x0) / sigma) ** 2)

def compton_energy(theta):
    return 511 / (2 - np.cos(theta))

# Create histogram of detected energies
counts, bin_edges = np.histogram(energies, bins=100)
bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

# Initial guess for parameters
a_init = max(counts)
x0_init = compton_energy(angle_rad)  
sigma_init = 0.0695 * compton_energy(angle_rad) / 2.355

# Create a LeastSquares object
cost = LeastSquares(bin_centers, counts, np.sqrt(counts), gaussian)

# Create a Minuit object with named parameters matching your function
m = Minuit(cost, a=a_init, x0=x0_init, sigma=sigma_init)

# Set parameter limits using the proper attribute syntax
m.limits["x0"] = (0, 511)
m.limits["sigma"] = (0, 511)

m.migrad()  # Perform the fit
# Extract fitted parameters
a_fit, x0_fit, sigma_fit = m.values[:3]
# Extract covariance matrix
covariance_matrix = m.covariance[:3, :3]
# Calculate the errors
errors = np.sqrt(np.diag(covariance_matrix))
# Calculate the correlation matrix
correlation_matrix = covariance_matrix / np.outer(errors, errors)
# Print the results
print("Fitted parameters:")
print(f"a: {a_fit} ± {errors[0]}")
print(f"x0: {x0_fit} ± {errors[1]}")
print(f"σ: {sigma_fit} ± {errors[2]}")
print("Correlation matrix:")
print(correlation_matrix)

# Plot the histogram and the fitted Gaussian
plt.figure(figsize=(12, 8))
plt.hist(energies, bins=100, alpha=0.4, color='blue', 
         edgecolor='black', label='Detected Energies')
plt.plot(bin_centers, gaussian(bin_centers, a_fit, x0_fit, sigma_fit),
            color='red', label='Fitted Gaussian')
plt.axvline(x0_init, color='green', linestyle='--', label='Compton Energy')

# Calculate additional statistics
FWHM = 2.355 * sigma_fit  # Full Width at Half Maximum
FWHM_err = 2.355 * errors[2]
ER = FWHM / x0_fit  # Energy Resolution
ER_err = np.sqrt((FWHM_err/x0_fit)**2 + (FWHM*errors[1]/x0_fit**2)**2)

# Calculate chi-squared
residuals = counts - gaussian(bin_centers, a_fit, x0_fit, sigma_fit)
chi2 = np.sum((residuals / np.sqrt(counts))**2)
ndf = len(counts) - 3  # number of data points minus number of parameters
chi2_ndf = chi2 / ndf

# Estimate number of events in peak (integral of Gaussian)
N_hit = a_fit * sigma_fit * np.sqrt(2*np.pi)
N_hit_err = np.sqrt((errors[0]*sigma_fit*np.sqrt(2*np.pi))**2 + 
                    (a_fit*errors[2]*np.sqrt(2*np.pi))**2)

# I = S * t * epsilon * solid_angle/4pi
# t = I / (S * epsilon * solid_angle/4pi)
S = 175000 * 903/1000 # Bq (Only for 511 KeV)
epsilon_gate = 0.4796 # Gate efficiency
solid_angle = 0.0197 # rad
n_run = len(file_names)  # Number of runs
I = N_tot * n_run  # beam intensity
time = 2 * I / (S * epsilon_gate * solid_angle / (4 * np.pi))  # s
# The factor 2 is because 511 keV is back to back

# Add statistics text box to the plot
stats_text = rf"$\chi^2/\mathrm{{ndf}} = {chi2:.1f}/{ndf} = {chi2_ndf:.2f}$" + "\n"
stats_text += f"<E> = {x0_fit:.2f} ± {errors[1]:.2f} keV\n"
stats_text += f"FWHM = {FWHM:.2f} ± {FWHM_err:.2f} keV\n"
stats_text += f"ER = {ER:.3f} ± {ER_err:.3f}\n"
stats_text += f"N = {N_hit:.1f} ± {N_hit_err:.1f}\n"
stats_text += f"Rate = {N_hit/time:.4f} ± {N_hit_err/time:.4f} Hz"

# Position the text box in the upper right corner
plt.text(0.3, 0.7, stats_text, transform=plt.gca().transAxes, fontsize=14, color="black", ha='left')

plt.title('Fitting Gaussian to Detected Energies')
plt.xlabel('Energy (keV)')
plt.ylabel('Counts')
plt.legend(fontsize=12, loc='upper left')
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.savefig(file_path + "plots/fitted_gaussian.png")

print(time)