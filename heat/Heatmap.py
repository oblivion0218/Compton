import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# Carica il file usando sep per spazi
title = "Calibrazione_Franco_AMP"
df = pd.read_csv(title +".txt", sep='\s+', header=None)

# Definire i limiti accettabili
lower_limit =  6.5 # Limite inferiore
upper_limit =  10 # Limite superiore

# Impostare i nomi delle colonne
df.columns = ["RE(%)", "sigma_RE(%)", "HV", "Gain"]

# Rimuovere la prima riga se contiene i nomi delle colonne come dati
df = df[1:]

# Converti le colonne numeriche a float, sostituendo i valori non numerici con NaN
df["RE(%)"] = pd.to_numeric(df["RE(%)"], errors='coerce')
df["HV"] = pd.to_numeric(df["HV"], errors='coerce')
df["Gain"] = pd.to_numeric(df["Gain"], errors='coerce')

# Rimuovere eventuali righe con NaN
df = df.dropna(subset=["RE(%)", "HV", "Gain"])

# Creazione della heatmap
heatmap_data = df.pivot(index="Gain", columns="HV", values="RE(%)")

# Creare una maschera per i valori accettabili e fuori scala
masked_data = heatmap_data.where((heatmap_data >= lower_limit) & (heatmap_data <= upper_limit))  # Valori validi
invalid_data = heatmap_data.where((heatmap_data < lower_limit) | (heatmap_data > upper_limit))  # Valori fuori scala

# Calcolare i limiti cromatici basati sui valori accettabili
vmin = masked_data.min().min()
vmax = masked_data.max().max()

# Configurazione e visualizzazione della heatmap
plt.figure(figsize=(10, 8))
ax = sns.heatmap(
    masked_data, annot=False, fmt=".2f", cmap="plasma",
    cbar_kws={'label': "Energetic resolution = "r"$\frac{FWHM}{<E_{\gamma}>} \cdot 100$"},
    linecolor='black', linewidth=0.5, vmin=vmin, vmax=vmax
)

# Annotazioni
for i in range(heatmap_data.shape[0]):
    for j in range(heatmap_data.shape[1]):
        value = heatmap_data.iloc[i, j]
        if pd.notnull(value):
            if lower_limit <= value <= upper_limit:
                # Annotazioni per i valori accettabili (bianco)
                ax.text(j + 0.5, i + 0.5, f"{value:.2f}", ha='center', va='center', color='white', fontsize=8)
            else:
                # Annotazioni per i valori fuori scala (nero)
                ax.text(j + 0.5, i + 0.5, f"{value:.2f}", ha='center', va='center', color='black', fontsize=8)

# Applicare lo sfondo bianco per i valori fuori scala
sns.heatmap(
    invalid_data.notnull(), cmap=["white"], cbar=False, linecolor='black', linewidths=0.5, ax=ax, mask=invalid_data.isnull()
)

ax.tick_params(axis='x', colors='black')  # Colore tick asse X
ax.tick_params(axis='y', colors='black')  # Colore tick asse Y
ax.spines['top'].set_color('black')  # Colore bordo superiore
ax.spines['bottom'].set_color('black')  # Colore bordo inferiore
ax.spines['left'].set_color('black')  # Colore bordo sinistro
ax.spines['right'].set_color('black')  # Colore bordo destro
plt.xlabel("HV (V)", color="black")
plt.ylabel("Gain", color="black")
plt.title(title, color="black")

# Titoli e assi
plt.xlabel("HV (V)")
plt.ylabel("Gain")
plt.title(title)
plt.savefig("heatmap_" + title + ".png", dpi=300)
plt.show()

