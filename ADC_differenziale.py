import ROOT
import MoraPyRoot as mpr
import matplotlib.pyplot as plt

file_path = "/mnt/c/Users/User/Desktop/info/Compton/ADC/"
# formatura = "gaussiana"
# formatura = "preamplificata"
formatura = "gaussiana-shaping"

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

channel_values, channel_errors = read_data(file_path + "dati/formatura_" + formatura + "_canale.txt")
counts_values, counts_errors = read_data(file_path + "dati/formatura_" + formatura + "_conteggi.txt")

plt.errorbar(channel_values, counts_values, xerr=channel_errors, yerr=counts_errors)

plt.xlabel("Canale")
plt.ylabel("Conteggi")
plt.ylim(0, 90)
plt.title("ADC differenziale formatura " + formatura)
plt.grid(True)
plt.savefig(file_path + "ADC_differenziale_formatura_" + formatura + ".png")
