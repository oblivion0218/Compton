import ROOT
import numpy as np
import math
from lib import MoraPyRoot as mpr
from lib import LabLibrary as ll

file_path = "/mnt/c/Users/User/Desktop/info/Compton/Measurments_trasmission/35_deg/"
file_path_background = "/mnt/c/Users/User/Desktop/info/Compton/Measurments_trasmission/50_deg/background/"
coo = [0.1, 0.63, 0.45, 0.9]

def create_hist(file_path, fileNamePNG):
    hs = ll.hist_vector(file_path + "data/")

    for h in hs:
        mpr.plot_hist_MPL(h, file_path + "plots/hist/" + h.GetName() + ".png")
    
    H = ll.spectum_sum(hs)
    mpr.plot_hist_MPL(H, file_path + "plots/hist/" + fileNamePNG)

    return H

#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
# I Approach - Fit peaks with background
#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
H = create_hist(file_path, "hist_sum.png")
integral_H = H.Integral()
H = ll.normalize_histogram(H)

def fit_peaks(hist, peak, extreme, fileNamePNG, graph_name, x_axis_name, y_axis_name, graphic_option, pave_coordinates = None):
    """
    Fit a peak with a gaussian function and a background with a linear function
    
    :param hist: ROOT histogram object.
    :param peak: Peak position.
    :param extreme: Range of the peak.
    :param fileNamePNG: Name of the output file.
    :param graph_name: Name of the graph.
    :param x_axis_name: Name of the x-axis. 
    :param y_axis_name: Name of the y-axis.
    :param graphic_option: Graphic option.
    :param pave_coordinates: Coordinates of the text box.
    """
    coo2 = [0.1, 0.5, 0.4, 0.9]
    str2 = ["f1", "f2", "Amp", "<x>", "#sigma"]

    n_bins = H.GetNbinsX()

    f_picco = ROOT.TF1("picco", "gaus(0)", 0, 2500)
    # L'inizializzazione dei parametri varia al variare del numero di bin
    f_picco.SetParameters(0, 40000 / n_bins)
    f_picco.SetParameters(1, peak)
    f_picco.SetParameters(2, 0.001 * n_bins)
    mpr.stampa_graph_fit(H, f_picco, file_path + "plots/fit/peak_" + fileNamePNG, "", "", "Counts", "", peak - extreme, peak + extreme)

    f_fondo = ROOT.TF1("fondo", "[0] + [1]/x", 0, 2500)
    f_fondo.SetParameters(0, 1)
    f_fondo.SetParameters(1, 1)
    mpr.stampa_graph_fit(H, f_fondo, file_path + "plots/fit/background_" + fileNamePNG, "", "", "Counts", "", peak - extreme * 3, peak + extreme * 3)

    # Background + peak
    f_true = ROOT.TF1("modello", "gaus(0) + [3] + [4]/x", 0, 2500)
    f_true.SetParameter(0, f_picco.GetParameter(0))
    f_true.SetParameter(1, f_picco.GetParameter(1))
    f_true.SetParameter(2, f_picco.GetParameter(2))
    f_true.SetParameter(3, f_fondo.GetParameter(0))
    f_true.SetParameter(4, f_fondo.GetParameter(1))

    min_val = peak - extreme
    max_val = peak + extreme

    ll.stampa_graph_fit_ComptonStudy(hist, f_true, integral_H, min_val, max_val, file_path, fileNamePNG, graph_name, 
                                     x_axis_name, y_axis_name, graphic_option, pave_coordinates, f_fondo)

peakCompton = ll.search_photopeak(H, 0.4, 4)
fit_peaks(H, peakCompton, 150, "Compton_fit.png", "Compton peak", "Energy [keV]", "Counts", "", coo)

peak511 = ll.search_photopeak(H, 0.2, 4)
fit_peaks(H, peak511, 100, "511_fit.png", "511 peak", "Energy [keV]", "Counts", "", coo)

#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
# II Approach - Fit peaks without background
#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
H_background = ll.normalize_histogram(create_hist(file_path_background, "hist_sum_background.png"))

# The idea is to normalize the background using the time of acquisition
# scale_factor = 9 / 18 # 50 deg
scale_factor = 9 / 14 # 35 deg

H_background = ll.normalize_histogram(H_background, scale_factor)

H_nb = ROOT.TH1F("hist_no_background", "Spectra without background", H.GetNbinsX(), H.GetXaxis().GetXmin(), H.GetXaxis().GetXmax())

for i in range(1, H.GetNbinsX() + 1):  # ROOT bins start from 1
    content = H.GetBinContent(i) - H_background.GetBinContent(i)
    H_nb.SetBinContent(i, content)
mpr.plot_hist_MPL(H_nb, file_path + "plots/hist/" + "hist_no_background.png")

def fit_peaks_no_background(hist, peak, extreme, fileNamePNG, graph_name, x_axis_name, y_axis_name, graphic_option, pave_coordinates = None):
    """
    Fit a peak with a gaussian function and a background with a linear function
    
    :param hist: ROOT histogram object.
    :param peak: Peak position.
    :param extreme: Range of the peak.
    :param fileNamePNG: Name of the output file.
    :param graph_name: Name of the graph.
    :param x_axis_name: Name of the x-axis. 
    :param y_axis_name: Name of the y-axis.
    :param graphic_option: Graphic option.
    :param pave_coordinates: Coordinates of the text box.
    """
    n_bins = H.GetNbinsX()

    f_true = ROOT.TF1("picco", "gaus(0)", 0, 2500)
    # L'inizializzazione dei parametri varia al variare del numero di bin
    f_true.SetParameters(0, 40000 / n_bins)
    f_true.SetParameters(1, peak)
    f_true.SetParameters(2, 0.001 * n_bins)

    min_val = peak - extreme
    max_val = peak + extreme

    ll.stampa_graph_fit_ComptonStudy(hist, f_true, integral_H, min_val, max_val, file_path, fileNamePNG, 
                                     graph_name, x_axis_name, y_axis_name, graphic_option, pave_coordinates)
    
peakCompton_nb = ll.search_photopeak(H_nb, 0.47, 4)
fit_peaks_no_background(H_nb, peakCompton_nb, 150, "Compton_fit_no_background.png", "Compton peak", "Energy [keV]", "Counts", "", coo)