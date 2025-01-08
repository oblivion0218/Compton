import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

file_path = "/mnt/c/Users/User/Desktop/info/Compton/Christmas_measurments/"

data = pd.read_csv(file_path + "results/Misure_invernali.txt", delim_whitespace=True)

# Estrae le colonne in array numpy
times_po = data["time"]
RE = data["RE"]
RE_error = data["sigma_RE"]
E_mean = data["<E>"]
E_mean_error = data["sigma_<E>"]

times = []
sum = 0
for t in times_po:
    t_str = str(t)
    n_hour = float(t_str[:1])  
    n_iteration = float(t_str[3:])  

    minutes = 60 * n_hour
    sum += minutes

    times.append(sum)

# Plot RE vs times con barre d'errore, punti più piccoli e griglia più fitta lungo l'asse x
plt.figure(figsize=(22, 5))
plt.errorbar(times, RE, yerr=RE_error, fmt='o', capsize=5, color='blue', markersize=2, label="Valori delle misure natalizie")  # Punti più piccoli
plt.xlabel('Time (min)')
plt.ylabel(r"$\frac{FWHM}{<E_{\gamma}>}$")
plt.title(r"$\frac{FWHM}{<E_{\gamma}>}$" + " vs Time")
plt.grid(True) 

true_value = 0.06952372784942827
true_error = 0.0003986994307091226
plt.axhline(true_value, 0, 23500, color="red", label="Valore misurato durante la calibrazione")
plt.fill_between(times, true_value - true_error, true_value + true_error, color='red', alpha=0.3)

plt.xticks(range(0, int(max(times)) + 10, 1000))
plt.legend(loc="lower right")

plt.savefig(file_path + "plots/RE_vs_time.png")

# Plot <E>_mean vs times con barre d'errore, punti più piccoli e griglia più fitta lungo l'asse x
plt.figure(figsize=(22, 5))
plt.errorbar(times, E_mean, yerr=E_mean_error, fmt='o', capsize=5, color='blue', markersize=2)  # Punti più piccoli
plt.xlabel('Time (min)')
plt.ylabel('<E>')
plt.title('<E> vs Time')
plt.grid(True)

plt.xticks(range(0, int(max(times)) + 10, 1000))

plt.savefig(file_path + "plots/<E>_vs_time.png")
