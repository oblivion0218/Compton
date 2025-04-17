import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.odr import ODR, Model, RealData

# --- Caricamento dati ---
file_path = ''
df = pd.read_csv(file_path + "Dati_calibrazione14_02_25.csv", delimiter=',', dtype=str)

df["ENERGIA"] = df["ENERGIA"].str.replace(',', '.').astype(float)
df["Canale associato"] = df["Canale associato"].str.replace(',', '.').astype(float)
df["err canale"] = df["err canale"].str.replace(',', '.').astype(float)

energies = df["ENERGIA"].values
channels = df["Canale associato"].values
errors = df["err canale"].values
elements = df["ELEMENTO"].values

# --- Modello lineare ---
def linear_model(B, x):
    return B[0] * x + B[1]

model = Model(linear_model)

# --- Dati con errore su x (canali) ---
data = RealData(channels, energies, sx=errors)  # sx = errore su x

# --- Fit con ODR ---
odr = ODR(data, model, beta0=[1., 0.])  # Stima iniziale: a=1, b=0
output = odr.run()

a, b = output.beta
a_err, b_err = output.sd_beta

# --- Plot ---
plt.figure(figsize=(8, 6))
plt.errorbar(channels, energies, xerr=errors, fmt='o', capsize=3, label='Data')

x_fit = np.linspace(min(channels), max(channels), 500)
y_fit = linear_model([a, b], x_fit)
plt.plot(x_fit, y_fit, color='red', label=fr'$y = {a:.2f}x + {b:.2f}$')

# Annotazioni
for i, txt in enumerate(elements):
    el = str(elements[i])[:2]
    n = str(elements[i])[-3:]
    label = fr"$^{{{n}}}{el}$"
    if i % 2 == 0:
        plt.annotate(label, (channels[i] + 50, energies[i] - 20), fontsize=12, ha='left')
    else:
        plt.annotate(label, (channels[i] - 30, energies[i]), fontsize=12, ha='right')

plt.xlabel("Channel")
plt.ylabel("Energy (keV)")
#plt.title("Energy vs. Channel with ODR Fit (Error on X)")
plt.grid()
plt.legend()
plt.savefig(file_path + "calibration14_02_25.png")

# --- Stampa dei risultati ---
print(f"Slope (a): {a:.5f} ± {a_err:.5f}")
print(f"Intercept (b): {b:.5f} ± {b_err:.5f}")
