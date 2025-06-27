import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# Carica il file usando sep per spazi
df = pd.read_csv("h_Franco_AMP.txt", sep='\s+', header=None)
title = "Franco - AMP"

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

# Imposta stile chiaro
sns.set_style("white")

# Crea figura e assi
fig, ax = plt.subplots(figsize=(10, 8))

# Heatmap con testo bianco all'interno
sns.heatmap(
    heatmap_data,
    annot=True,
    fmt=".1f",
    cmap="plasma",
    linecolor='black',
    linewidth=0.5,
    annot_kws={"color": "white", "weight": "bold", "size": 12},
    ax=ax
)

# Etichette assi e titolo in nero
ax.set_xlabel("HV (V)", fontsize=14, fontweight='bold', color='black')
ax.set_ylabel("Gain", fontsize=14, fontweight='bold', color='black')

# Tick assi in nero
ax.tick_params(axis='both', colors='black', labelsize=12)

# Colorbar in nero
cbar = ax.collections[0].colorbar
cbar.ax.yaxis.label.set_color('black')
cbar.ax.tick_params(colors='black', labelsize=12)

# Layout finale e salvataggio
plt.tight_layout()
plt.savefig("PP_" + title + ".png", dpi=300, bbox_inches='tight')
plt.show()

