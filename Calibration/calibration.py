import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from iminuit import Minuit
from iminuit.cost import LeastSquares

# Load the CSV file (make sure to replace 'your_file.csv' with the actual filename)
file_path = '/mnt/c/Users/User/Desktop/info/Compton/Calibration/'
df = pd.read_csv(file_path + "Dati_calibrazione - Foglio1.csv", delimiter=',', dtype=str)

# Convert numeric columns, replacing commas with dots
df["ENERGIA"] = df["ENERGIA"].str.replace(',', '.').astype(float)
df["Canale associato"] = df["Canale associato"].str.replace(',', '.').astype(float)
df["err canale"] = df["err canale"].str.replace(',', '.').astype(float)

# Extract data
energies = df["ENERGIA"].values
channels = df["Canale associato"].values
errors = df["err canale"].values
elements = df["ELEMENTO"].values

# Define linear model
def linear_model(x, a, b):
    return a * x + b

# Define least squares cost function
least_squares = LeastSquares(channels, energies, errors, linear_model)
m = Minuit(least_squares, a=1, b=0)
m.migrad()  # Run the fit

# Plot data points with error bars
plt.figure(figsize=(8, 6))
plt.errorbar(channels, energies, yerr=errors, fmt='o', label='Data', capsize=3)

# Plot best fit line
x_fit = np.linspace(min(channels), max(channels), 500)
y_fit = linear_model(x_fit, *m.values)
plt.plot(x_fit, y_fit, label=f'$y = {m.values[0]:.2f}x {m.values[1]:.2f}$', color='red')

# Annotate points
for i, txt in enumerate(elements):
    el = str(elements[i])[:2]
    n = str(elements[i])[-3:]

    label = fr"$^{{{n}}}{el}$"

    if i % 2 == 0:
        plt.annotate(label, (channels[i] + 50, energies[i] - 20), fontsize=12, ha='left')
    else:
        plt.annotate(label, (channels[i] - 30, energies[i]), fontsize=12, ha='right')

# Labels and title
plt.xlabel("Channel")
plt.ylabel("Energy (keV)")
plt.title("Energy vs. Channel with Linear Fit")
plt.legend()
plt.grid()
plt.savefig(file_path + "calibration.png")

# Print fit results
print(m)
