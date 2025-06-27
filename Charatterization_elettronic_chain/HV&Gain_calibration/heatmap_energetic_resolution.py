import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# Carica il file
title = "_Franco_AMP-TSCA"
df = pd.read_csv("calibrazione" + title + ".txt", sep='\s+', header=None)

# Limiti accettabili
lower_limit = 6.5
upper_limit = 10

# Nomi colonne
df.columns = ["RE(%)", "sigma_RE(%)", "HV", "Gain"]
df = df[1:]

# Conversioni
df["RE(%)"] = pd.to_numeric(df["RE(%)"], errors='coerce')
df["HV"] = pd.to_numeric(df["HV"], errors='coerce')
df["Gain"] = pd.to_numeric(df["Gain"], errors='coerce')
df = df.dropna(subset=["RE(%)", "HV", "Gain"])

# Pivot per heatmap
heatmap_data = df.pivot(index="Gain", columns="HV", values="RE(%)")
masked_data = heatmap_data.where((heatmap_data >= lower_limit) & (heatmap_data <= upper_limit))
invalid_data = heatmap_data.where((heatmap_data < lower_limit) | (heatmap_data > upper_limit))
vmin = masked_data.min().min()
vmax = masked_data.max().max()

# Plot
plt.figure(figsize=(10, 8))
ax = sns.heatmap(
    masked_data, annot=False, fmt=".2f", cmap="plasma",
    linecolor='black', linewidth=0.5, vmin=vmin, vmax=vmax
)

# Annotazioni: ora testo bianco/nero + **GRASSETTO**
for i in range(heatmap_data.shape[0]):
    for j in range(heatmap_data.shape[1]):
        value = heatmap_data.iloc[i, j]
        if pd.notnull(value):
            is_valid = lower_limit <= value <= upper_limit
            ax.text(j + 0.5, i + 0.5, f"{value:.2f}",
                    ha='center', va='center',
                    color='white' if is_valid else 'black',
                    fontsize=12, fontweight='bold')  # GRASSETTO QUI

# Sfondo bianco dove i valori sono fuori range
sns.heatmap(
    invalid_data.notnull(), cmap=sns.color_palette(["white"], as_cmap=True),
    cbar=False, linecolor='black', linewidths=0.5, ax=ax,
    mask=invalid_data.isnull()
)

# Ticks e bordi in nero
ax.tick_params(axis='x', colors='black')
ax.tick_params(axis='y', colors='black')
for spine in ax.spines.values():
    spine.set_color('black')

# Titoli assi con **GRASSETTO**
ax.set_xlabel("HV (V)", color="black", fontsize=14, fontweight='bold')
ax.set_ylabel("Gain", color="black", fontsize=14, fontweight='bold')

# Salvataggio
plt.tight_layout()
plt.savefig("ER" + title + ".png", dpi=300)
plt.show()
