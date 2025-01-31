import ROOT
import numpy as np
import MoraPyRoot as mpr
import LabLibrary as ll

#Andrea
file_path = "/mnt/c/Users/User/Desktop/info/Compton/Christmas_measurments/"

# ritorna il valore della risoluzione energetica (FWHM / <E>) e l'energia media al picco (<E>)
def fit_photopeak(hist, fileNamePNG, noise_threshold, n_peaks):
    
    f_true = ll.fit_photopeak_linear_background(hist, fileNamePNG, noise_threshold, n_peaks)[1]

    FWHM = 2.3548 * f_true.GetParameter(4)
    mean = f_true.GetParameter(3)
    mean_err = f_true.GetParError(3) 
    risoluzione_energetica =  FWHM / mean
    sigma_risoluzione_energetica = np.sqrt((2.3548 * f_true.GetParError(4) / f_true.GetParameter(3)) ** 2 + (FWHM *  f_true.GetParError(3) / (f_true.GetParameter(3) ** 2)) ** 2)
    


    return ((risoluzione_energetica, sigma_risoluzione_energetica), (mean, mean_err))

#-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

fileNames = ll.remove_extension(ll.get_file_names(file_path + "data_gennaio_4/"))

RE = []
E_mean = []
times = []

for fileName in fileNames:
    data = ll.read_histogram_data(file_path + "data_gennaio_4/" + fileName + ".Spe")

    hist = ROOT.TH1D("h", "h", len(data), 0, len(data))
    time = fileName[7:]

    for i in range(len(data)):
        hist.SetBinContent(i + 1, data[i])  # i+1 perché i bin in ROOT partono da 1

    mpr.plot_hist_MPL(hist, file_path + "plots/hist/h_gennaio_4_" + time + ".png")
    fit_result = fit_photopeak(hist, file_path + "plots/fit/h_fit_gennaio_4_" + time + ".png", noise_threshold=0.5, n_peaks=2)

    RE.append(fit_result[0])
    E_mean.append(fit_result[1])
    times.append(time)

with open("Misure_gennaio_4.txt", "w") as file:
    file.write("\nRE\tsigma_RE\t<E>\tsigma_<E>\ttime\n")
    for r, e, t in zip(RE, E_mean, times):
        file.write(f"{r[0]}\t{r[1]}\t{e[0]}\t{e[1]}\t{t}\n")

print("\nFine\n")
