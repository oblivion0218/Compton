import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.optimize import curve_fit
from scipy import stats
from iminuit import Minuit
from iminuit.cost import LeastSquares

# Add a new constant at the top of your script:
N_tot = 1000000  # Total number of photons in each simulation

# Set the directory containing the simulation result files
data_dir = "../simulated_events/data"

# Initialize data structures to store metrics across simulation runs
sim_runs = []
photons_left_target = []
interacting_photons = []
photoelectric_interactions = []
compton_1st = []
compton_2nd = []
compton_3rd = []
photons_reached_detector = []
photons_observed = []
all_energies = []

# Get all the txt files in the directory without using glob
files = []
for file_name in os.listdir(data_dir):
    if file_name.endswith("_detected_energies.txt"):
        files.append(os.path.join(data_dir, file_name))
files.sort()

# Process each file
for file_path in tqdm(files, desc="Processing files"):
    # Extract run number from filename
    run_number = int(os.path.basename(file_path).split("_")[0])
    sim_runs.append(run_number)
    
    # Parse the metrics from the file
    with open(file_path, 'r') as f:
        lines = f.readlines()
        
        # Extract metrics
        photons_left_target.append(int(lines[0].split(": ")[1]))
        interacting_photons.append(int(lines[1].split(": ")[1]))
        photoelectric_interactions.append(int(lines[2].split(": ")[1]))
        compton_1st.append(int(lines[3].split(": ")[1]))
        compton_2nd.append(int(lines[4].split(": ")[1]))
        compton_3rd.append(int(lines[5].split(": ")[1]))
        photons_reached_detector.append(int(lines[6].split(": ")[1]))
        photons_observed.append(int(lines[7].split(": ")[1]))
        
        # Extract energy values, starting from line 9
        energies = [float(line.strip()) for line in lines[9:] if line.strip()]
        all_energies.extend(energies)

# Convert lists to numpy arrays for easier handling
sim_runs = np.array(sim_runs)
photons_left_target = np.array(photons_left_target)
interacting_photons = np.array(interacting_photons)
photoelectric_interactions = np.array(photoelectric_interactions)
compton_1st = np.array(compton_1st)
compton_2nd = np.array(compton_2nd)
compton_3rd = np.array(compton_3rd)
photons_reached_detector = np.array(photons_reached_detector)
photons_observed = np.array(photons_observed)
all_energies = np.array(all_energies)


plt.figure(figsize=(12, 8))

# Plotting all metrics on the same plot
plt.plot(sim_runs, interacting_photons/N_tot, '-o', color='red', label='Photon Interactions')
plt.plot(sim_runs, photoelectric_interactions/N_tot, '-o', color='green', label='Photoelectric Interactions')
plt.plot(sim_runs, compton_1st/N_tot, '-o', color='purple', label='Compton 1st')
plt.plot(sim_runs, compton_2nd/N_tot, '-o', color='orange', label='Compton 2nd')
plt.plot(sim_runs, compton_3rd/N_tot, '-o', color='brown', label='Compton 3rd')

plt.title('Photon Interaction probability across simulation runs')
plt.xlabel('Simulation Run')
plt.ylabel('Probability / Scaled Efficiency')
plt.legend(fontsize=12)
plt.grid(True)
plt.savefig("interaction_probabilities.png")

plt.figure(figsize=(12, 8))
plt.plot(sim_runs, photons_reached_detector/N_tot, '-o', color='cyan', label='Photons Reaching Detector')
plt.plot(sim_runs, photons_observed/N_tot, '-o', color='magenta', label='Photons Observed')

detection_efficiency = 100 * photons_observed / photons_reached_detector

sorted_indices = np.argsort(sim_runs)
sim_runs_sorted = sim_runs[sorted_indices]
detection_efficiency_sorted = detection_efficiency[sorted_indices]

plt.plot(sim_runs, detection_efficiency / 100, '-o', color='black', label='Detection Efficiency (scaled)')

plt.title('Photon interaction across first simulation runs')
plt.xlabel('Simulation run')
plt.ylabel('Probability / Scaled Efficiency')
plt.legend(fontsize=12)
plt.grid(True)
plt.savefig("detector_parameter.png")

# Histogram of detected photon energies
plt.figure(figsize=(12, 8))
plt.hist(all_energies, bins=100, alpha=0.7, color='blue', edgecolor='black', density=True)
plt.title('Probability Density of Detected Photon Energies')
plt.xlabel('Energy (keV)')
plt.ylabel('Probability Density')
plt.grid(True)

# Calculate statistics for the energy distribution
mean_energy = np.mean(all_energies)
median_energy = np.median(all_energies)
std_energy = np.std(all_energies)

# Add statistical information to the plot
plt.annotate(f'Mean: {mean_energy:.2f} keV\nMedian: {median_energy:.2f} keV\nStd Dev: {std_energy:.2f} keV',
             xy=(0.7, 0.85), xycoords='axes fraction',
             bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="black", lw=1))

# Define Gaussian function for fitting
def gaussian(x, amplitude, mean, sigma):
    return amplitude * np.exp(-(x - mean)**2 / (2 * sigma**2))

# Create histogram of all detected energies as probability density
plt.figure(figsize=(12, 8))

# First plot the histogram
counts, bin_edges, _ = plt.hist(all_energies, bins=100, alpha=0.7, color='blue', 
                               edgecolor='black', density=True)
bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])

# Initial guess for the Gaussian parameters
# Assume the photopeak is around the median energy
max_count_idx = np.argmax(counts)
initial_amplitude = counts[max_count_idx]
initial_mean = bin_centers[max_count_idx]
initial_sigma = std_energy / 2  # Start with a reasonable guess

# Define a region around the peak to fit (to avoid fitting the Compton continuum)
# This assumes the photopeak is the highest point in the histogram
peak_region_mask = (bin_centers > initial_mean - 3*initial_sigma) & (bin_centers < initial_mean + 3*initial_sigma)
x_peak = bin_centers[peak_region_mask]
y_peak = counts[peak_region_mask]

try:
    # Fit the Gaussian to the peak region
    popt, pcov = curve_fit(gaussian, x_peak, y_peak, 
                          p0=[initial_amplitude, initial_mean, initial_sigma])
    
    # Extract the parameters
    amplitude, mean_energy_fit, sigma = popt
    
    # Calculate FWHM (Full Width at Half Maximum)
    fwhm = 2.355 * sigma  # FWHM = 2.355 * sigma for a Gaussian
    
    # Calculate energy resolution
    energy_resolution = fwhm / mean_energy_fit

    #calculating the rate
    rate = 0
    
    # Count photons in the photopeak (within ±2 sigma)
    photopeak_mask = (all_energies > mean_energy_fit - 2*sigma) & (all_energies < mean_energy_fit + 2*sigma)
    n_photopeak = np.sum(photopeak_mask)
    
    # Calculate percentage of photons in the photopeak
    percentage_in_photopeak = 100 * n_photopeak / len(all_energies)
    
    # Generate x values for the fitted curve
    x_fit = np.linspace(bin_edges[0], bin_edges[-1], 1000)
    y_fit = gaussian(x_fit, *popt)
    
    # Plot the fitted curve
    plt.plot(x_fit, y_fit, 'r-', linewidth=2, 
            label='Gaussian fit\n'
                  f'Mean: {mean_energy_fit:.2f} keV\n'
                  f'FWHM: {fwhm:.2f} keV\n'
                  f'Resolution: {energy_resolution*100:.2f}%\n'
                  f'Photons in peak: {n_photopeak}\n'
                  f'Peak fraction: {percentage_in_photopeak:.2f}%'
                  f'Peak Rate: {rate:.2f}%')
    
    # Highlight the FWHM on the plot
    half_max = amplitude / 2.0
    left_idx = np.abs(y_fit - half_max).argmin()
    while x_fit[left_idx] > mean_energy_fit and left_idx > 0:
        left_idx -= 1
    while y_fit[left_idx] < half_max and left_idx > 0:
        left_idx -= 1
    
    right_idx = np.abs(y_fit - half_max).argmin()
    while x_fit[right_idx] < mean_energy_fit and right_idx < len(x_fit) - 1:
        right_idx += 1
    while y_fit[right_idx] < half_max and right_idx < len(x_fit) - 1:
        right_idx += 1
    
    plt.plot([x_fit[left_idx], x_fit[right_idx]], [half_max, half_max], 'g-', linewidth=2, 
            label='FWHM')
    
    # Add vertical lines at the FWHM points
    plt.axvline(x=x_fit[left_idx], color='g', linestyle='--', alpha=0.5)
    plt.axvline(x=x_fit[right_idx], color='g', linestyle='--', alpha=0.5)
    
    plt.legend(loc='upper right')
    
except RuntimeError:
    print("Error: Could not fit Gaussian to the photopeak")
    # Add fallback plot with just the statistics we already have
    plt.annotate(f'Mean: {mean_energy:.2f} keV\nMedian: {median_energy:.2f} keV\nStd Dev: {std_energy:.2f} keV',
                xy=(0.7, 0.85), xycoords='axes fraction',
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="black", lw=1))

plt.title('Energy Distribution with Gaussian Fit of Photopeak')
plt.xlabel('Energy (keV)')
plt.ylabel('Probability Density')
plt.grid(True)

# Save the updated histogram with Gaussian fit
plt.tight_layout()
plt.savefig(os.path.join(data_dir, 'new_energy_histogram_with_fit.png'))

# Create histogram of all detected energies
plt.figure(figsize=(12, 8))

# First plot the histogram
counts, bin_edges, _ = plt.hist(all_energies, bins=100, alpha=0.7, color='blue', 
                               edgecolor='black', density=True)
bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])

# Define Gaussian function for fitting
def gaussian(x, amplitude, mean, sigma):
    return amplitude * np.exp(-(x - mean)**2 / (2 * sigma**2))

# For 511 keV source scattered at 50 degrees, expected energy is approximately:
expected_energy = 511 / (1 + (511/511) * (1 - np.cos(np.radians(50))))
print(f"Expected energy for 511 keV photon scattered at 50°: {expected_energy:.2f} keV")

# Initial guess: Look for a peak near the expected energy
nearby_mask = np.abs(bin_centers - expected_energy) < 50  # Look within ±50 keV
if np.sum(nearby_mask) > 0:
    local_max_idx = np.argmax(counts[nearby_mask])
    max_count_idx = np.where(nearby_mask)[0][local_max_idx]
else:
    # Fall back to global maximum if no peak found near expected energy
    max_count_idx = np.argmax(counts)

initial_amplitude = counts[max_count_idx]
initial_mean = bin_centers[max_count_idx]
initial_sigma = 20  # Start with a reasonable guess for a detector

# Define a region around the peak for fitting
peak_region_mask = (bin_centers > initial_mean - 4*initial_sigma) & (bin_centers < initial_mean + 4*initial_sigma)
x_peak = bin_centers[peak_region_mask]
y_peak = counts[peak_region_mask]

if len(x_peak) > 5:  # Need enough points for a meaningful fit
    try:
        # Define the least squares cost function
        least_squares = LeastSquares(x_peak, y_peak, np.sqrt(y_peak), gaussian)
        
        # Initialize Minuit
        m = Minuit(least_squares, amplitude=initial_amplitude, mean=initial_mean, sigma=initial_sigma)
        
        # Set limits (optional, but helpful)
        m.limits["amplitude"] = (0, None)
        m.limits["sigma"] = (0, None)
        m.limits["mean"] = (initial_mean - 100, initial_mean + 100)
        
        # Perform the fit
        m.migrad()
        m.hesse()  # Get more accurate error estimates
        
        # Extract the parameters
        amplitude = m.values["amplitude"]
        mean_energy_fit = m.values["mean"]
        sigma = m.values["sigma"]
        
        # Calculate FWHM (Full Width at Half Maximum)
        fwhm = 2.355 * sigma  # FWHM = 2.355 * sigma for a Gaussian
        
        # Calculate energy resolution
        energy_resolution = fwhm / mean_energy_fit
        
        # Count photons in the photopeak (within ±2 sigma)
        photopeak_mask = (all_energies > mean_energy_fit - 2*sigma) & (all_energies < mean_energy_fit + 2*sigma)
        n_photopeak = np.sum(photopeak_mask)
        
        # Calculate percentage of photons in the photopeak
        percentage_in_photopeak = 100 * n_photopeak / len(all_energies)
        
        # Generate x values for the fitted curve
        x_fit = np.linspace(bin_edges[0], bin_edges[-1], 1000)
        y_fit = gaussian(x_fit, amplitude, mean_energy_fit, sigma)
        
        # Plot the fitted curve
        plt.plot(x_fit, y_fit, 'r-', linewidth=2, 
                label='Gaussian fit\n'
                    f'Mean: {mean_energy_fit:.2f} keV\n'
                    f'FWHM: {fwhm:.2f} keV\n'
                    f'Resolution: {energy_resolution*100:.2f}%\n'
                    f'Photons in peak: {n_photopeak}\n'
                    f'Peak fraction: {percentage_in_photopeak:.2f}%')
        
        # Highlight the FWHM on the plot
        half_max = amplitude / 2.0
        
        # Find left and right points where the curve crosses half-maximum
        left_x = mean_energy_fit - sigma * np.sqrt(2 * np.log(2))
        right_x = mean_energy_fit + sigma * np.sqrt(2 * np.log(2))
        
        plt.plot([left_x, right_x], [half_max, half_max], 'g-', linewidth=2, label='FWHM')
        
        # Add vertical lines at the FWHM points
        plt.axvline(x=left_x, color='g', linestyle='--', alpha=0.5)
        plt.axvline(x=right_x, color='g', linestyle='--', alpha=0.5)
        
        # Add a vertical line at the theoretical value
        plt.axvline(x=expected_energy, color='purple', linestyle='-', alpha=0.5,
                   label=f'Expected: {expected_energy:.2f} keV')
        
        plt.legend(loc='upper right', fontsize=12)
        
        # Add detailed fit results to console output
        print("\nPhotopeak Analysis (iminuit):")
        print("-" * 40)
        print(f"Fit converged: {m.valid}")
        print(f"Fit quality: {m.fmin}")
        print(f"Expected energy (511 keV @ 50°): {expected_energy:.2f} keV")
        print(f"Fitted mean energy: {mean_energy_fit:.2f} ± {m.errors['mean']:.2f} keV")
        print(f"Fitted sigma: {sigma:.2f} ± {m.errors['sigma']:.2f} keV")
        print(f"FWHM: {fwhm:.2f} keV")
        print(f"Energy resolution (FWHM/E): {energy_resolution*100:.2f}%")
        print(f"Number of photons in photopeak: {n_photopeak}")
        print(f"Percentage of detected photons in photopeak: {percentage_in_photopeak:.2f}%")
        print(f"Deviation from expected: {(mean_energy_fit - expected_energy)/expected_energy*100:.2f}%")
        
    except Exception as e:
        print(f"iminuit fitting error: {e}")
        plt.title('Energy Distribution (Fitting Failed)')
else:
    print("Not enough data points near the expected energy peak")
    plt.title('Energy Distribution (Insufficient Data)')

plt.xlabel('Energy (keV)')
plt.ylabel('Probability Density')
plt.grid(True)

# Save the updated histogram with Gaussian fit
plt.tight_layout()
plt.savefig(os.path.join(data_dir, 'energy_histogram_with_fit.png'))
