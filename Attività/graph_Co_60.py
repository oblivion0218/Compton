import numpy as np
import matplotlib.pyplot as plt

# Funzione per leggere il file
def read_data(filename):
    data = []
    with open(filename, 'r') as f:
        for line in f:
            if line.strip() and not line.startswith("Canale"):  # Ignora intestazione
                try:
                    parts = line.split()
                    if len(parts) >= 3:  # Controllo che ci siano almeno 3 valori
                        _, _, eventi = map(int, parts[-3:])  # Prende l'ultimo valore come conteggio
                        data.append(eventi)
                except ValueError:
                    continue  # Ignora righe non valide
    return np.array(data)

# Lettura dei dati
filename = "misura_Co60_SOLOsorgente.txt"
data = read_data(filename)

# Creazione dell'istogramma
plt.figure(figsize=(10, 5))
plt.hist(np.arange(len(data)), bins=len(data), weights=data, color='b', alpha=0.7, label="Conteggi")
plt.xlabel("Canale")
plt.ylabel("Conteggi")
plt.legend()
plt.grid()
plt.show()
