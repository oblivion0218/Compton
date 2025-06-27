import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

data = pd.read_csv('fit_results.txt', delim_whitespace=True, skiprows=1, 
                   names=["x", "Offset", "Offset_err", "Slope", "Slope_err", 
                          "Amplitude", "Amplitude_err", "Center", "Center_err", 
                          "Sigma", "Sigma_err"])

# Calcola la risoluzione e la propagazione degli errori
data['Resolution'] = (2.355 *100* data['Sigma']) / data['Center']
# Propagazione dell'errore per R = (2.48 * Sigma) / Center
data['Resolution_err'] = data['Resolution'] * np.sqrt(
    (data['Sigma_err'] / data['Sigma'])**2 + (data['Center_err'] / data['Center'])**2)

# Grafico
plt.errorbar(data['x'], data['Resolution'], yerr=data['Resolution_err'], fmt='o', capsize=3, label="Risoluzione")
plt.plot(data['x'], data['Resolution'], linestyle='-', color='C0')  # unisce i punti con una linea
plt.xlabel("Shaping time (ns)")
plt.ylabel("Resolution (%)")
plt.grid()
plt.savefig("Risoluzione-shaping-FRANCO.png")