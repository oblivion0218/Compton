import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# Carica il file usando sep per spazi
df = pd.read_csv("h_Ugo_2_AMP.txt", sep='\s+', header=None)
title = "Ugo - AMP"

# Impostare i nomi delle colonne
df.columns = ["h_peak", "err_h", "HV", "Gain"]

# Rimuovere la prima riga se contiene i nomi delle colonne come dati
df = df[1:]

# Converti le colonne numeriche a float, sostituendo i valori non numerici con NaN
df["h_peak"] = pd.to_numeric(df["h_peak"], errors='coerce')
df["HV"] = pd.to_numeric(df["HV"], errors='coerce')
df["Gain"] = pd.to_numeric(df["Gain"], errors='coerce')

# Rimuovere eventuali righe con NaN
df = df.dropna(subset=["h_peak", "HV", "Gain"])

# Creazione della heatmap
heatmap_data = df.pivot(index="Gain", columns="HV", values="h_peak")

# Configurazione e visualizzazione della heatmap
plt.figure(figsize=(10, 8))
sns.heatmap(
    heatmap_data, annot=heatmap_data, fmt=".1f", cmap="plasma",
    cbar_kws={'label': "Peak position"},
    linecolor='black', linewidth=0.5, annot_kws={"color": "black"}
)
plt.xlabel("HV (V)")
plt.ylabel("Gain")
plt.title(title)
plt.savefig("h_" + title + ".png", dpi=300)