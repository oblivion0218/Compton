import ROOT
import numpy as np
import MoraPyRoot as mpr
import LabLibrary as ll

#Andrea
file_path = "/mnt/c/Users/User/Desktop/info/Compton/HV&Gain_calibration/AMP/CalibrazioneRivelatoreUgo/"

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
def fit_photopeak(hist, fileNamePNG, noise_threshold, n_peaks):
    
    f_true = ll.fit_photopeak_linear_background(hist, fileNamePNG, noise_threshold, n_peaks)[1]

    FWHM = 2.3548 * f_true.GetParameter(4)
    risoluzione_energetica =  FWHM / f_true.GetParameter(3)
    sigma_risoluzione_energetica = np.sqrt((2.3548 * f_true.GetParError(4) / f_true.GetParameter(3)) ** 2 + (FWHM *  f_true.GetParError(3) / (f_true.GetParameter(3) ** 2)) ** 2)
    
    return (risoluzione_energetica, sigma_risoluzione_energetica)

#-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

fileNames = ll.remove_extension(ll.get_file_names(file_path + "HV&Gain/"))

RE = []

for fileName in fileNames:
    voltage = int(fileName.split('V_')[0])
    gain = int(fileName.split('V_')[1].split('G')[0])

    data_brutti = ll.read_histogram_data(file_path + "HV&Gain/" + fileName + ".Spe")
    data = pulizia_dati(data_brutti, MaxNZeros=100, noise_threshold=50) # [Ugo (50, 50) Franco (30, 100)](AMP) 

    # if fileName == "550V_500G_0F": # correzione per AMP-TSCA Ugo
    #     data = pulizia_dati(data_brutti, MaxNZeros=400, noise_threshold=50) # [Ugo (50, 50) Franco (30, 100)](AMP) 
    # else: 
    #     data = pulizia_dati(data_brutti, MaxNZeros=100, noise_threshold=50) # [Ugo (50, 50) Franco (30, 100)](AMP) 

    # Creazione istogramma
    hist = ROOT.TH1D("h", "h", len(data), 0, len(data))

    for i in range(len(data)):
        hist.SetBinContent(i + 1, data[i])  # i+1 perché i bin in ROOT partono da 1

    mpr.plot_hist_MPL(hist, file_path + fileName + ".png")
    RE.append((fit_photopeak(hist, file_path + "h_fit_" + fileName + ".png", noise_threshold=0.5, n_peaks=2), voltage, gain))

with open("Calibrazione_Ugo_AMP.txt", "w") as file:
    file.write("\nRE(%)\tsigma_RE(%)\tHV\tGain\n")
    for r in RE:
        file.write(f"{r[0][0] * 100}\t{r[0][1] * 100}\t{r[1]}\t{r[2]}\n")

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

# # mpr.plot_hist_MPL(hist, file_path + fileName + ".png")
# print(fit_photopeak(hist, file_path + fileName + ".png", noise_threshold=0.3, n_peaks=2))

#--------------------------------------------------------------------------------------------
# Franco - AMP-TSCA (600V, 10G) (550V, 20G) (500V, 50G)
#--------------------------------------------------------------------------------------------

# RE_bis = []

# fileNames = ["600V_10G_0F", "550V_20G_0F", "500V_50G_0F"]

# for fileName in fileNames:
#     voltage = int(fileName.split('V_')[0])
#     gain = int(fileName.split('V_')[1].split('G')[0])

#     data_brutti = ll.read_histogram_data( file_path + "HV&Gain/" + fileName + ".Spe")
#     data = pulizia_dati(data_brutti, MaxNZeros=100, noise_threshold=50) # [Ugo (50, 50) Franco (30, 100)](AMP) 

#     # Creazione istogramma
#     hist = ROOT.TH1D("h", "h", 1500, 0, 1500)

#     for i in range(1500):
#         hist.SetBinContent(i + 1, data[i])  # i+1 perché i bin in ROOT partono da 1

#     RE_bis.append((fit_photopeak(hist, file_path + fileName + ".png", noise_threshold=0.6, n_peaks=2), voltage, gain))
# print(RE_bis)

#--------------------------------------------------------------------------------------------
# Ugo - AMP-TSCA (750V, 20G) (850V, 10G)
#--------------------------------------------------------------------------------------------

# RE_bis = []

# fileNames = ["850V_10G_0F", "750V_20G_0F"]

# for fileName in fileNames:
#     voltage = int(fileName.split('V_')[0])
#     gain = int(fileName.split('V_')[1].split('G')[0])

#     data_brutti = ll.read_histogram_data( file_path + "HV&Gain/" + fileName + ".Spe")
#     data = pulizia_dati(data_brutti, MaxNZeros=100, noise_threshold=50) # [Ugo (50, 50) Franco (30, 100)](AMP) 

#     # Creazione istogramma
#     hist = ROOT.TH1D("h", "h", 1500, 0, 1500)

#     for i in range(1500):
#         hist.SetBinContent(i + 1, data[i])  # i+1 perché i bin in ROOT partono da 1

#     RE_bis.append((fit_photopeak(hist, file_path + fileName + ".png", noise_threshold=0.6, n_peaks=2), voltage, gain))
# print(RE_bis)
