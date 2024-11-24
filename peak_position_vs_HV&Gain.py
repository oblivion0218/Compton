import ROOT
import numpy as np
import MoraPyRoot as mpr
import LabLibrary as ll

#Andrea
percorso_file = "/mnt/c/Users/User/Desktop/info/AMP-TSCA/CalibrazioneRivelatoreUgo/"

# Funzione per cancellare i canali finali con solo 0
def pulizia_dati(data, MaxNZeros, noise_threshold=0):
    x = []
    r = 0

    # Rimuove i primi e gli ultimi canali rumorosi
    for i in range(len(data)):
        if i <= 20 or i > 2000:
            x.append(0)
        else:
            x.append(data[i])

        # Conta zeri consecutivi
        if data[i] <= noise_threshold:
            r += 1
            if i > 0 and data[i - 1] != 0:
                r = 0
        if r > MaxNZeros:
            break

    # # Rimuove zeri iniziali e crea il risultato finale
    # x_def = []
    # start = False
    # for val in x:
    #     if val != 0:
    #         start = True
    #     if start:
    #         x_def.append(val)
            
    return x

# ritorna il valore della risoluzione energetica = FWHM / <E>
def fit_photopeak(hist, fileName, noise_threshold, n_peaks):
    # Definizione delle coordinate e delle etichette per le legende
    coo1 = [0.1, 0.6, 0.4, 0.9]  # coordinate legenda
    str1 = ["Amp", "<x>", "#sigma"]  # etichette legenda
    coo2 = [0.6, 0.13, 0.9, 0.53]
    str2 = ["f1", "f2", "Amp", "<x>", "#sigma"]

    photopeak_x = ll.search_photopeak(hist, noise_threshold, n_peaks)
    
    n_bins = hist.GetNbinsX()

    # Determina l'estremo per il fit
    if n_bins >= 1800: 
        extreme = hist.GetXaxis().GetBinUpEdge(n_bins) / 10
    else: 
        extreme = hist.GetXaxis().GetBinUpEdge(n_bins) / 50

    # Fit gaussiano
    f_picco = ROOT.TF1("picco", "gaus(0)", 0, 2000)
    # L'inizializzazione dei parametri varia al variare del numero di bin
    f_picco.SetParameters(0, 40000 / n_bins)
    f_picco.SetParameters(1, photopeak_x)
    f_picco.SetParameters(2, 0.001 * n_bins)
    mpr.fit(hist, f_picco, 3, "S", 1000, photopeak_x - extreme, photopeak_x + extreme)

    # Fit lineare
    f_fondo = ROOT.TF1("fondo", "pol1(0)", 0, 2000)
    f_fondo.SetParameters(0, 1)
    f_fondo.SetParameters(1, 1)
    mpr.fit(hist, f_fondo, 3, "S", 1000, photopeak_x - extreme * 2, photopeak_x + extreme * 2)

    # Modello combinato (fondo + picco)
    f_true = ROOT.TF1("modello", "pol1(0) + gaus(2)", 0, 2000)
    f_true.SetParameter(0, f_fondo.GetParameter(0))
    f_true.SetParameter(1, f_fondo.GetParameter(1))
    f_true.SetParameter(2, f_picco.GetParameter(0))
    f_true.SetParameter(3, f_picco.GetParameter(1))
    f_true.SetParameter(4, f_picco.GetParameter(2))
    mpr.stampa_graph_fit(hist, f_true, photopeak_x - extreme, photopeak_x + extreme, percorso_file + "h_fit" + fileName + ".png", "", "", "Counts", "", photopeak_x - extreme, photopeak_x + extreme, 5)

    return (f_true.GetParameter(3),  f_true.GetParError(3))

#-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

fileNames = ll.remove_extension(ll.get_file_names(percorso_file + "HV&Gain/"))

RE = []

for fileName in fileNames:
    voltage = int(fileName.split('V_')[0])
    gain = int(fileName.split('V_')[1].split('G')[0])

    data_brutti = ll.read_histogram_data(percorso_file + "HV&Gain/" + fileName + ".Spe")
    data = pulizia_dati(data_brutti, MaxNZeros=100, noise_threshold=50) # [Ugo (50, 50) Franco (30, 100)](AMP) 

    # if fileName == "550V_500G_0F": # correzione per AMP-TSCA Ugo
    #     data = pulizia_dati(data_brutti, MaxNZeros=400, noise_threshold=50) # [Ugo (50, 50) Franco (30, 100)](AMP) 
    # else: 
    #     data = pulizia_dati(data_brutti, MaxNZeros=100, noise_threshold=50) # [Ugo (50, 50) Franco (30, 100)](AMP) 

    # Creazione istogramma
    hist = ROOT.TH1D("h", "h", len(data), 0, len(data))

    for i in range(len(data)):
        hist.SetBinContent(i + 1, data[i])  # i+1 perché i bin in ROOT partono da 1

    mpr.plot_hist_MPL(hist, percorso_file + fileName + ".png")
    RE.append((fit_photopeak(hist, fileName, noise_threshold=0.5, n_peaks=2), voltage, gain))

with open("h_Ugo_AMP-TSCA.txt", "w") as file:
    file.write("\nh_peak\terr_h\tHV\tGain\n")
    for r in RE:
        file.write(f"{r[0][0]}\t{r[0][1]}\t{r[1]}\t{r[2]}\n")

print("\nFine\n")

#-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
# Refit canali problematico 
#-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
# Franco - AMP (550V, 100G)
#--------------------------------------------------------------------------------------------

# data_brutti = ll.read_histogram_data(percorso_file + "HV&Gain/550V_100G_0F.Spe")
# data = pulizia_dati(data_brutti, MaxNZeros=30, noise_threshold=100) # Ugo (50, 50) Franco (30, 100)

# # Creazione istogramma
# hist = ROOT.TH1D("h", "h", len(data), 0, len(data))
# for i in range(len(data)):
#     hist.SetBinContent(i + 1, data[i])  # i+1 perché i bin in ROOT partono da 1

# # mpr.plot_hist_MPL(hist, percorso_file + fileName + ".png")
# print(fit_photopeak(hist, "550V_100G_0F", noise_threshold=0.3, n_peaks=2))

#--------------------------------------------------------------------------------------------
# Franco - AMP-TSCA (600V, 10G) (550V, 20G) (500V, 50G)
#--------------------------------------------------------------------------------------------

# RE_bis = []

# fileNames = ["600V_10G_0F", "550V_20G_0F", "500V_50G_0F"]

# for fileName in fileNames:
#     voltage = int(fileName.split('V_')[0])
#     gain = int(fileName.split('V_')[1].split('G')[0])

#     data_brutti = ll.read_histogram_data( percorso_file + "HV&Gain/" + fileName + ".Spe")
#     data = pulizia_dati(data_brutti, MaxNZeros=100, noise_threshold=50) # [Ugo (50, 50) Franco (30, 100)](AMP) 

#     # Creazione istogramma
#     hist = ROOT.TH1D("h", "h", 1500, 0, 1500)

#     for i in range(1500):
#         hist.SetBinContent(i + 1, data[i])  # i+1 perché i bin in ROOT partono da 1

#     RE_bis.append((fit_photopeak(hist, fileName, noise_threshold=0.6, n_peaks=2), voltage, gain))
# print(RE_bis)

#--------------------------------------------------------------------------------------------
# Ugo - AMP-TSCA (750V, 20G) (850V, 10G)
#--------------------------------------------------------------------------------------------

# RE_bis = []

# fileNames = ["850V_10G_0F", "750V_20G_0F"]

# for fileName in fileNames:
#     voltage = int(fileName.split('V_')[0])
#     gain = int(fileName.split('V_')[1].split('G')[0])

#     data_brutti = ll.read_histogram_data( percorso_file + "HV&Gain/" + fileName + ".Spe")
#     data = pulizia_dati(data_brutti, MaxNZeros=100, noise_threshold=50) # [Ugo (50, 50) Franco (30, 100)](AMP) 

#     # Creazione istogramma
#     hist = ROOT.TH1D("h", "h", 1500, 0, 1500)

#     for i in range(1500):
#         hist.SetBinContent(i + 1, data[i])  # i+1 perché i bin in ROOT partono da 1

#     RE_bis.append((fit_photopeak(hist, fileName, noise_threshold=0.6, n_peaks=2), voltage, gain))
# print(RE_bis)
