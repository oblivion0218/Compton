import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

file_path = "/mnt/c/Users/User/Desktop/info/Compton/Christmas_measurments/"

# Read data from the text file using whitespace as delimiter
data = pd.read_csv(file_path + "results/Misure_gennaio_3.txt", delim_whitespace=True)

# Extract specific columns into numpy arrays
times_po = data["time"]
RE = data["RE"]
RE_error = data["sigma_RE"]
E_mean = data["<E>"]
E_mean_error = data["sigma_<E>"]

# Convert time to minutes and accumulate the time values
times = []
sum = 0
for t in times_po:
    t_str = str(t)

    # Misure invernali

    # n_hour = float(t_str[:1])  # Extract hour value
    # n_iteration = float(t_str[3:])  # Extract iteration part of the time

    # minutes = 60 * n_hour  # Convert hours to minutes
    # sum += minutes  # Accumulate the total time in minutes

    # # Misure gennaio
    # minutes = float(t_str[:4]) / 60 # Extract hour value 20 min
    # n_iteration = float(t_str[8:])  # Extract iteration part of the time

    # Misure gennaio 2
    minutes = 20

    sum += minutes
    times.append(sum)

# Plot RE vs times
plt.figure(figsize=(22, 5))
plt.errorbar(times, RE, yerr=RE_error, fmt='o', capsize=5, color='blue', markersize=2, label="Values from Christmas measurements")  # Small points for data
plt.xlabel('Time (min)')
plt.ylabel(r"$\frac{FWHM}{<E_{\gamma}>}$")  # Label for FWHM over gamma energy
plt.title(r"$\frac{FWHM}{<E_{\gamma}>}$" + " vs Time")  # Plot title
plt.grid(True)  # Enable grid for the plot

# Add a horizontal line for the measured calibration value
true_value = 0.06952372784942827
true_error = 0.0003986994307091226
plt.axhline(true_value, 0, 23500, color="red", label="Measured value during calibration")
plt.fill_between(times, true_value - true_error, true_value + true_error, color='red', alpha=0.3)

# Set x-ticks and plot the legend
plt.xticks(range(0, int(max(times)) + 10, 1000))
plt.legend(loc="lower right")
plt.savefig(file_path + "plots/RE_vs_time_gennaio_3.png")  # Save the figure to a file

# Plot <E>_mean vs times 
plt.figure(figsize=(22, 5))
plt.errorbar(times, E_mean, yerr=E_mean_error, fmt='o', capsize=5, color='blue', markersize=2 , label="Values from Christmas measurements")  # Small points for data
plt.xlabel('Time (min)')
plt.ylabel('<E>')  # Energy mean
plt.title('<E> vs Time')  # Plot title
plt.grid(True)

# Add a horizontal line for the measured calibration value
# true_value = 1398.5585559377448
# true_error = 0.27309155059643514
# plt.axhline(true_value, 0, 23500, color="red", label="Measured value during calibration")
# plt.fill_between(times, true_value - true_error, true_value + true_error, color='red', alpha=0.3)

# Set x-ticks and plot the legend
plt.xticks(range(0, int(max(times)) + 10, 1000))
plt.legend(loc="lower left")
plt.savefig(file_path + "plots/<E>_vs_time_gennaio_3.png")  # Save the figure to a file

# #*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*- 

# # Load temperature data
# data_temp = pd.read_csv(file_path + "temperature_data.csv", sep=";")

# # Clean column names by stripping spaces
# data_temp.columns = data_temp.columns.str.strip()

# # Create list of time in minutes
# minuti = [i for i in range(len(data_temp["Ora"]))]
# temperature = data_temp["Temperatura_Celsius"]

# # Create the plots
# fig, axs = plt.subplots(3, 2, figsize=(22, 15), sharex=True, gridspec_kw={'height_ratios': [3, 2]})

# # First plot: <E>_mean vs time with error bars
# axs[0].errorbar(times, E_mean, yerr=E_mean_error, fmt='o', capsize=5, color='blue', markersize=2, label="Values from Christmas measurements")
# axs[0].set_ylabel('<E>')
# axs[0].set_title('<E> vs Time')  # Plot title
# axs[0].grid(True)
# axs[0].legend(loc="lower left")

# # Second plot: temperature over time
# axs[1].plot(minuti, temperature, color='red', label='Temperature (°C)')
# axs[1].set_xlabel('Time (min)')
# axs[1].set_ylabel('Temperature (°C)')  # Temperature in Celsius
# axs[1].grid(True)
# axs[1].legend()

# # Set x-ticks for the time axis in minutes
# plt.xticks(range(0, int(max(minuti)) + 10, 1000))

# # Save the plot to a file
# plt.savefig(file_path + "plots/<E>_vs_time_with_temperature.png")
