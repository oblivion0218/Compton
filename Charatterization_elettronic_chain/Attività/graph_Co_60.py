import numpy as np
import matplotlib.pyplot as plt

# Funzione per leggere il file
import numpy as np

def read_data_fixed(filename):
    data = []
    with open(filename, 'r', encoding="utf-8") as f:
        reading_data = False  # Flag per iniziare la lettura solo dopo l'intestazione
        for line in f:
            if line.startswith("Canale"):  # Riconosce l'inizio dei dati
                reading_data = True
                continue
            if reading_data:
                parts = line.strip().split("\t")  # Usa tabulazione come separatore
                if len(parts) == 3:  # Assicura di avere tre colonne
                    try:
                        eventi = int(parts[2])  # Prende l'ultima colonna (Eventi)
                        data.append(eventi)
                    except ValueError:
                        continue  # Ignora eventuali righe non valide
    return np.array(data)

# Lettura dei dati
filename = "misura_Co60_SOLOsorgente.txt"# Provo a leggere i dati corretti
data_fixed = read_data_fixed(filename)

energia_min, energia_max = 1100, 1400

# Creazione dell'istogramma
plt.figure(figsize=(10, 5))
plt.hist(np.arange(len(data_fixed)), bins=len(data_fixed), weights=data_fixed, color='b', alpha=0.7, label="Counts")
plt.xlabel("Channel")
plt.ylabel("Counts")
plt.legend()
plt.grid()

min_valore , max_valore = 600 , 800
plt.xlim(min_valore, max_valore)


plt.show()
plt.savefig("Doppietto.png")
