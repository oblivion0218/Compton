import numpy as np
import matplotlib.pyplot as plt

# Inizializzazione liste
angle = []
rate = []
err_rate = []
N = []
err_N = []
centr = []
err_centr = []

file_path = "/mnt/c/Users/ASUS/Desktop/WSL_shared/Compton/Codes/data_analysis/"

# Lettura dati dal file (saltando la prima riga)
with open(file_path + "parameters.txt", "r") as file:
    next(file)  # salta intestazione
    for line in file:
        values = line.strip().split()
        if len(values) == 7:
            angle.append(float(values[0]))
            rate.append(float(values[1]))
            err_rate.append(float(values[2]))
            N.append(float(values[3]))
            err_N.append(float(values[4]))
            centr.append(float(values[5]))
            err_centr.append(float(values[6]))

# --- Ordinamento in base all'angolo ---
combined = list(zip(angle, rate, err_rate, N, err_N, centr, err_centr))
combined.sort(key=lambda x: x[0])  # ordina per angolo (prima colonna)

# Riassegna i valori ordinati
angle, rate, err_rate, N, err_N, centr, err_centr = zip(*combined)

# Funzione per creare e salvare i plot
def crea_plot(x, y, yerr, ylabel, filename):
    plt.figure()
    plt.errorbar(x, y, yerr=yerr, fmt='o-', capsize=4, label=ylabel)
    plt.xlabel("Angle (deg)")
    plt.ylabel(ylabel)
    plt.title(f"{ylabel} vs Angle")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()

# Creazione e salvataggio dei plot
crea_plot(angle, rate, err_rate, "Rate", "plot_rate.png")
crea_plot(angle, N, err_N, "Counts", "plot_Counts.png")
crea_plot(angle, centr, err_centr, "Center", "plot_center.png")

print("Plot salvati --> DONE ")
