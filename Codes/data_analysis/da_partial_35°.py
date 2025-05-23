import ROOT
import numpy as np
import matplotlib.pyplot as plt
from lib import MoraPyRoot as mpr
from lib import LabLibrary as ll


file_path = "/mnt/c/Users/User/Desktop/info/Compton/Measurments/Measurments_trasmission/35_deg/"
# file_path = "/mnt/c/Users/ASUS/Desktop/WSL_shared/Compton/Measurments/Measurments_trasmission/35_deg/" 


def fit_peaks(hist, peak, sigma, left_step, right_step, x_axis_name, y_axis_name, file_path):
    """
    Fit a peak with a gaussian function and a background with a linear function 
    
    :param hist: ROOT histogram object.
    :param peak: Peak position.
    :param sigma: Range of the peak.
    :param left_step: Left step from the center of the Compton peak for the fit.
    :param right_step: Right step from the center of the Compton peak for the fit.
    :param x_axis_name: Name of the x-axis. 
    :param y_axis_name: Name of the y-axis.
    :param file_path: Path to save the plots.
    """
    coo0 = [0.1, 0.65, 0.45, 0.9]
    str0 = ["Amp", "<x>", "#sigma"]
    coo1 = [0.1, 0.5, 0.45, 0.9]
    coo2 = [0.1, 0.35, 0.45, 0.9]
    
    #-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    # PARTIAL FIT - Compton peak
    #-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    f_Compton = ROOT.TF1("Compton_peak", "gaus(0)", 0, 2000)
    f_Compton.SetParameter(0, 100) 
    f_Compton.SetParameter(1, peak)
    f_Compton.SetParameter(2, sigma)

    mpr.stampa_graph_fit(hist, f_Compton, file_path + "Compton_peak_.png", "Compton peak", 
                         x_axis_name, y_axis_name, "", peak - 2 * sigma, peak + 2 * sigma, 3, coo0, str0)
    
    #-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    # PARTIAL FIT - 511 peak
    #-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    f_511 = ROOT.TF1("511_peak", "gaus(0)", 0, 2000)
    f_511.SetParameter(0, 100) 
    f_511.SetParameter(1, 1500)
    f_511.SetParameter(2, sigma/2)

    mpr.stampa_graph_fit(hist, f_511, file_path + "511_peak_.png", "Compton peak", 
                         x_axis_name, y_axis_name, "", 1500 - 2 * sigma, 1500 + 2 * sigma, 3, coo0, str0)
    
    #-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    # PARTIAL FIT - background
    #-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    f_background = ROOT.TF1("background", "pol4(0) + gaus(5)", 0, 2000)
    f_background.SetParameter(0, 1)
    f_background.SetParameter(1, 1)
    f_background.SetParameter(2, 1)
    f_background.SetParameter(3, 1)
    f_background.SetParameter(4, 1)
    f_background.FixParameter(5, f_511.GetParameter(0))
    f_background.FixParameter(6, f_511.GetParameter(1))
    f_background.FixParameter(7, f_511.GetParameter(2))

    mpr.stampa_graph_fit(hist, f_background, file_path + "background_.png", "Background",
                         x_axis_name, y_axis_name, "", 600, 1800)

    #-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    # FIT - Complete model
    #-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    f_true = ROOT.TF1("model", "pol4(0) + gaus(5) + gaus(8)", 0, 2000) 
    f_true.SetParameter(0, f_background.GetParameter(0))
    f_true.SetParameter(1, f_background.GetParameter(1))
    f_true.SetParameter(2, f_background.GetParameter(2))
    f_true.SetParameter(3, f_background.GetParameter(3))
    f_true.SetParameter(4, f_background.GetParameter(4))
    f_true.SetParameter(5, f_Compton.GetParameter(0))
    f_true.SetParameter(6, f_Compton.GetParameter(1))
    f_true.SetParameter(7, f_Compton.GetParameter(2))
    f_true.SetParameter(8, f_background.GetParameter(5))
    f_true.SetParameter(9, f_background.GetParameter(6))
    f_true.SetParameter(10, f_background.GetParameter(7))

    fit_result = mpr.stampa_graph_fit(hist, f_true, file_path + "final_fit_.png", "Spectrum", x_axis_name, y_axis_name, 
                                      "", f_Compton.GetParameter(1) - left_step, f_Compton.GetParameter(1) + right_step)

    f_background.SetParameter(0, f_true.GetParameter(0))
    f_background.SetParameter(1, f_true.GetParameter(1))    
    f_background.SetParameter(2, f_true.GetParameter(2))
    f_background.SetParameter(3, f_true.GetParameter(3))
    f_background.SetParameter(4, f_true.GetParameter(4))
    f_background.SetParameter(5, f_true.GetParameter(8))
    f_background.SetParameter(6, f_true.GetParameter(9))
    f_background.SetParameter(7, f_true.GetParameter(10))

    return fit_result, f_background, f_true


#-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# Main 
#-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
angle = 35
time = 43000 * 14

H = ll.create_hist(file_path, "hist_sum.png")
peakCompton = ll.search_photopeak(H, 0.4, 2, file_path + "plots/fit/find_Compton_peak.png")
sigmaCompton = 50

left_step = 250
right_step = 200

hist_integral = H.Integral()
rebin_param = 3
H.Rebin(rebin_param)

fit_result, f_background, f_true = fit_peaks(H, peakCompton, sigmaCompton, left_step, right_step, "Energy [channels]", "Counts", 
                                             file_path + "plots/fit/")

E_mean = (f_true.GetParameter(6), f_true.GetParError(6))
sigma = (f_true.GetParameter(7), f_true.GetParError(7))
min_fit = E_mean[0] - left_step
max_fit = E_mean[0] + right_step

# Final fit
counts , rate = ll.plot_results(H, fit_result, f_background, f_true, rebin_param, time, E_mean, sigma, min_fit, max_fit, 
                                file_path + "plots/fit/", "fit_results.png", "Energy [channels]", "Counts")

ll.update_or_append_line("parameters.txt", angle, rate[0], rate[1], counts[0], counts[1], E_mean[0], E_mean[1], sigma[0], sigma[1])
