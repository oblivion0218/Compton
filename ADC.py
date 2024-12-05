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

amp_values, amp_errors = read_data(file_path + "dati/formatura_" + formatura + "_amp(mV).txt")
channel_values, channel_errors = read_data(file_path + "dati/formatura_" + formatura + "_canale.txt")
counts_values, counts_errors = read_data(file_path + "dati/formatura_" + formatura + "_conteggi.txt")

# Convert data in array for ROOT 
import array
amp_values_array = array.array('d', amp_values)
channel_values_array = array.array('d', channel_values)
amp_errors_array = array.array('d', amp_errors)
channel_errors_array = array.array('d', channel_errors)

graph_amp = ROOT.TGraphErrors(len(amp_values), amp_values_array, channel_values_array, amp_errors_array, channel_errors_array)

f_true = ROOT.TF1("ADC", "pol1(0)", 0, 2000)
f_true.SetParameters(0, 0)
f_true.SetParameters(1, 0.2)

coo1 = [0.1, 0.6, 0.45, 0.9]  
str1 = ["q", "m"] 
mpr.stampa_graph_fit(graph_amp, f_true, file_path +  "ADC_integrale_formatura_" + formatura + ".png", "ADC integrale formatura " + formatura, "Ampiezza (mV)", "Canale", "AP", 0, 7500, 2, coo1, str1)

plt.errorbar(channel_values, counts_values, xerr=channel_errors, yerr=counts_errors)

plt.xlabel("Canale")
plt.ylabel("Conteggi")
plt.ylim(0, 90)
plt.title("ADC differenziale formatura " + formatura)
plt.grid(True)
plt.savefig(file_path + "ADC_differenziale_formatura_" + formatura + ".png")
