import numpy as np
import matplotlib.pyplot as plt

# Percorso del file caricato
file_path = "CoincidenceON.Spe"

# Funzione sistemata per leggere i dati dal file
def read_data_fixed(filename):
    data = []
    with open(filename, 'r', encoding="utf-8") as f:
        reading_data = False
        for line in f:
            if line.startswith("$DATA"):
                reading_data = True
                continue
            if reading_data:
                try:
                    parts = line.strip()
                    eventi = int(parts)
                    data.append(eventi)
                except ValueError:
                    continue
    return np.array(data)

# Lettura dei dati
data = read_data_fixed(file_path)

# Creazione dell'istogramma
plt.figure(figsize=(10, 6))
#plt.title("Istogramma dei dati da CoincidenceON.Spe")
plt.hist(np.arange(len(data)), bins=len(data), weights=data, color='b', alpha=0.7, label="Counts")
plt.xlabel("Channel")
plt.ylabel("Counts")
plt.legend()

min_valore , max_valore = 400 , 600
plt.xlim(min_valore, max_valore)


plt.show()
plt.savefig("coincidenceON.png")
