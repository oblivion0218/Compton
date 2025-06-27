import ROOT
import matplotlib.pyplot as plt

file_path = "/mnt/c/Users/ASUS/Desktop/WSL_Shared/Compton/Charatterization_elettronic_chain/ADC/"
formatura = ["gaussiana", "preamplificata", "gaussiana-shaping"]
title = ["Gaussian", "Preamplified", "Gaussian + shaping"]

# Funzione per leggere i dati da un file
def read_data(fileName):
    x_values, x_errors = [], []
    with open(fileName, 'r') as file:
        for line in file:
            parts = line.strip().split()
            if len(parts) == 2:
                x, dx = map(float, parts)
                x_values.append(x)
                x_errors.append(dx)
    return x_values, x_errors

plt.figure(figsize=(16, 6))
i = 0
for name in formatura:
    channel_values, channel_errors = read_data(file_path + "data/formatura_" + name + "_canale.txt")
    counts_values, counts_errors = read_data(file_path + "data/formatura_" + name + "_conteggi.txt")
    plt.errorbar(channel_values, counts_values, xerr=channel_errors, yerr=counts_errors, label = title[i])
    i+= 1

plt.xlabel("Channel", fontsize=14)
plt.ylabel("Counts", fontsize=14)
plt.tick_params(labelsize=12)  # ingrandisce numeri sugli assi
plt.ylim(50, 90)
plt.grid(True)
plt.legend(fontsize=16)
plt.savefig(file_path + "plots/ADC_differenziale_TOT.png")
