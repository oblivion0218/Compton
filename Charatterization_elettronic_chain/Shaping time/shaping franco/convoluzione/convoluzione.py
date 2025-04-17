import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

data = pd.read_csv('convoluzione_shaping.txt', delim_whitespace=True, skiprows=1, 
                   names=["N","centro","errcentro",	"FWHM", "errFWHM"])

# Calcola la risoluzione e la propagazione degli errori
data['Resolution'] =  100*data['FWHM'] / data['centro']
data['Resolution_err'] = data['Resolution'] * np.sqrt(
    (data['errcentro'] / data['centro'])**2 + (data['errFWHM'] / data['FWHM'])**2)

# Grafico
plt.errorbar(data['N'], data['Resolution'], yerr=data['Resolution_err'], fmt='o', capsize=3, label="Risoluzione")
plt.xlabel("Shaping time (ns)")
plt.ylabel("Risoluzione")
plt.title("Risoluzione con CONVOLUZIONE")
plt.savefig("Risoluzione-shaping-FRANCO-CONVOLUZIONE.png")