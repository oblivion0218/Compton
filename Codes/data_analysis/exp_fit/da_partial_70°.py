import ROOT
import numpy as np
import matplotlib.pyplot as plt
from lib import MoraPyRoot as mpr
from lib import LabLibrary as ll


#file_path = "/mnt/c/Users/User/Desktop/info/Compton/Measurments/Measurments_riflection/70_deg_2nd/"
file_path = "/mnt/c/Users/ASUS/Desktop/WSL_shared/Compton/Measurments/Measurments_riflection/70_deg_2nd/"

def fit_peaks(hist, peak, sigma, min_fit, max_fit, x_axis_name, y_axis_name, file_path):
    """
    Fit a peak with a gaussian function and a background with a linear function 
    
    :param hist: ROOT histogram object.
    :param peak: Peak position.
    :param sigma: Range of the peak.
    :param min_fit: Minimum value of the fit.
    :param max_fit: Maximum value of the fit.
    :param x_axis_name: Name of the x-axis. 
    :param y_axis_name: Name of the y-axis.
    :param file_path: Path to save the plots.
    """
    coo0 = [0.55, 0.65, 0.9, 0.9]
    str0 = ["Amp", "<x>", "#sigma"]
    coo1 = [0.55, 0.5, 0.9, 0.9]
    coo2 = [0.55, 0.35, 0.9, 0.9]

    #-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    # PARTIAL FIT - Background
    #-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    f_back_e = ROOT.TF1("f_background", "expo(0)", 0, 2000)
    f_back_e.SetParameter(0, 1)
    f_back_e.SetParameter(1, 1)

    mpr.stampa_graph_fit(hist, f_back_e, file_path + "background_exp_.png", "Compton peak", 
                         x_axis_name, y_axis_name, "", 550, 1400, 2, coo0, ["f1", "f2"])
    
    #-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    # PARTIAL FIT - Compton peak
    #-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    f_Compton = ROOT.TF1("Compton_peak", "gaus(0)", 0, 2000)
    f_Compton.SetParameter(0, 100) 
    f_Compton.SetParameter(1, peak)
    f_Compton.SetParameter(2, sigma)

    mpr.stampa_graph_fit(hist, f_Compton, file_path + "Compton_peak_.png", "Compton peak", 
                         x_axis_name, y_axis_name, "", peak - 3 * sigma, peak + 3 * sigma, 3, coo0, str0)
    
    #-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    # FIT - Complete model
    #-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    
    f_true = ROOT.TF1("model", "expo(0) + gaus(2)", 0, 2000) 
    f_true.SetParameter(0, f_back_e.GetParameter(0))
    f_true.SetParameter(1, f_back_e.GetParameter(1))
    f_true.SetParameter(2, f_Compton.GetParameter(0))
    f_true.SetParameter(3, f_Compton.GetParameter(1))
    f_true.SetParameter(4, f_Compton.GetParameter(2))

    fit_result = mpr.stampa_graph_fit(hist, f_true, file_path + "final_fit_.png", "Spectrum", 
                                      x_axis_name, y_axis_name, "", min_fit, max_fit, 5, coo2, ["f1", "f2", "Amp", "<x>", "#sigma"])

    f_back_e.SetParameter(0, f_true.GetParameter(0))
    f_back_e.SetParameter(1, f_true.GetParameter(1))

    return fit_result, f_back_e, f_true


# #-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# # Main 
# #-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
H = ll.create_hist(file_path, "hist_sum.png")
peakCompton = ll.search_photopeak(H, 0.4, 2, file_path + "plots/fit/find_Compton_peak.png")
sigmaCompton = 50

# Study of the stability of the fit changing fit extremes, in order to choose the domain of the fit
step = 30
max_step = 20
# ll.stability_study_extreme(fit_peaks, H, peakCompton, sigmaCompton, step, max_step, "Energy [channels]", "Counts", file_path + "plots/fit/")

n_steps = 4
min_fit = peakCompton - (n_steps + 3) * step
max_fit = peakCompton + n_steps * step

# Study of the stability of the fit changing hist rebin, in order to choose the domain of the fit
rebin_max = 28
# ll.stability_study_rebin(fit_peaks, H, peakCompton, sigmaCompton, rebin_max, min_fit, max_fit, "Energy [channels]", "Counts", file_path + "plots/fit/")

hist_integral = H.Integral()
rebin_param = 3
H.Rebin(rebin_param)

fit_result, f_background, f_true = fit_peaks(H, peakCompton, sigmaCompton, min_fit, max_fit, "Energy [channels]", "Counts", 
                                             file_path + "plots/fit/")


time = 16 * 43000

# Final fit
counts , rate = ll.plot_results(H, fit_result, f_background, f_true, rebin_param, min_fit, max_fit, file_path + "plots/fit/", 
                "fit_results.png", "Energy [channels]", "Counts", time)

centroid = f_true.GetParameter(3)
centroid_err = f_true.GetParError(3)
angle = 70

ll.update_or_append_line("parameters.txt", angle, rate[0], rate[1], counts[0], counts[1], centroid, centroid_err)
