import pandas as pd 

NomeFile = "ADC - Foglio1"
tab = pd.read_csv(NomeFile + ".csv").drop(columns="CANALE 1","CANALE 2")

print(tab)