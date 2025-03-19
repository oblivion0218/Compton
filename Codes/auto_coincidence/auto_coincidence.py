from lib import MoraPyRoot as mpr
from lib import LabLibrary as ll
import ROOT

file_path = "/mnt/c/Users/User/Desktop/info/Compton/Auto-coincidence/"

def fit_photopeak(hist, fileNamePNG, noise_threshold, n_peaks):
    f_true = ll.fit_photopeak_linear_background(hist, fileNamePNG, noise_threshold, n_peaks)[1]
    return (f_true.GetParameter(3),  f_true.GetParError(3))

#-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

fileNames = ll.remove_extension(ll.get_file_names(file_path + "data"))
h = {}

for fileName in fileNames:

    data = ll.read_histogram_data(file_path + "data/" + fileName + ".Spe") 
    
    hist = ROOT.TH1D("h", "h", len(data), 0, len(data))

    for i in range(len(data)):
        hist.SetBinContent(i + 1, data[i])  # i+1 perché i bin in ROOT partono da 1

    mpr.plot_hist_MPL(hist, file_path + fileName + ".png")
    h[fileName] = fit_photopeak(hist, file_path + "h_fit_" + fileName + ".png", noise_threshold=0.5, n_peaks=2)

print(h)

