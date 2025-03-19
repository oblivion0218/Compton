import ROOT
import numpy as np
import matplotlib.pyplot as plt
from lib import MoraPyRoot as mpr
from lib import LabLibrary as ll


file_path = "/mnt/c/Users/User/Desktop/info/Compton/Measurments_trasmission/50_deg/"


def create_hist(file_path, fileNamePNG):
    hs = ll.hist_vector(file_path + "data/")

    for h in hs:
        mpr.plot_hist_MPL(h, file_path + "plots/hist/" + h.GetName() + ".png")
    
    H = ll.spectum_sum(hs)
    mpr.plot_hist_MPL(H, file_path + "plots/hist/" + fileNamePNG)

    return H


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
    coo0 = [0.1, 0.65, 0.45, 0.9]
    str0 = ["Amp", "<x>", "#sigma"]
    coo1 = [0.1, 0.5, 0.45, 0.9]
    coo2 = [0.1, 0.35, 0.45, 0.9]

    #-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    # PARTIAL FIT - Background
    #-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    f_back_e = ROOT.TF1("f_background", "expo(0)", 0, 2000)
    f_back_e.SetParameter(0, 1)
    f_back_e.SetParameter(1, 1)

    mpr.stampa_graph_fit(hist, f_back_e, file_path + "background_exp_.png", "Compton peak", 
                         x_axis_name, y_axis_name, "", 500, 2000, 2, coo0, ["f1", "f2"])
    
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


def stability_study_extreme(H, peakCompton, sigmaCompton, step, max_step, x_axis_name, y_axis_name, file_path):
    """"
    Study of the stability of the fit changing fit extremes, in order to choose the domain of the fit

    :param H: ROOT histogram object.
    :param peakCompton: Peak position.
    :param sigmaCompton: Range of the peak.
    :param step: Step of the fit extremes.
    :param max_step: Maximum number of steps.
    :param x_axis_name: Name of the x-axis.
    :param y_axis_name: Name of the y-axis.
    :param file_path: Path to save the plots.
    """
    chi2 = []
    n_steps = [i for i in range(1, max_step)]

    for i in n_steps:
        min_fit = peakCompton - i * step
        max_fit = peakCompton + i * step
        c,_,_= fit_peaks(H, peakCompton, sigmaCompton, min_fit, max_fit, x_axis_name, y_axis_name, 
                         file_path + "fit_stability/" + str(i))
        chi2.append(c.Chi2() / c.Ndf())

    plt.figure()
    plt.plot(n_steps, chi2, color="blue", marker="o", linestyle="dashed")
    plt.xlabel(r"$N_{step}$")
    plt.ylabel(r"$\chi^2 / ndf$")
    plt.grid()
    plt.xticks(n_steps)
    plt.savefig(file_path + "chi2_n_step.png")
    plt.close()


def stability_study_rebin(H, peakCompton, sigmaCompton, rebin_index, min_fit, max_fit, x_axis_name, y_axis_name, file_path):
    """
    Study of the stability of the fit changing hist rebin, in order to choose the domain of the fit

    :param H: ROOT histogram object.
    :param peakCompton: Peak position.
    :param sigmaCompton: Range of the peak.
    :param rebin_index: List of rebin factors.
    :param min_fit: Minimum value of the fit.
    :param max_fit: Maximum value of the fit.
    :param x_axis_name: Name of the x-axis.
    :param y_axis_name: Name of the y-axis.
    :param file_path: Path to save the
    """
    chi2 = []

    for r in rebin_index: 
        h = H.Clone()
        h.Rebin(r)
        h.Scale(1 / r)
        c, _, _ = fit_peaks(h, peakCompton, sigmaCompton, min_fit, max_fit, x_axis_name, y_axis_name, 
                            file_path + "fit_stability/Rebin_" + str(r))
        chi2.append(c.Chi2() / c.Ndf())

    plt.figure()
    plt.plot(rebin_index, chi2, color="blue", marker="o", linestyle="dashed")
    plt.xlabel(r"$Rebin factor$")
    plt.ylabel(r"$\chi^2 / ndf$")
    plt.grid()
    plt.xticks(rebin_index)
    plt.savefig(file_path + "chi2_rebin.png")
    plt.close()


#-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# Main 
#-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
H = create_hist(file_path, "hist_sum.png")
peakCompton = ll.search_photopeak(H, 0.4, 2, file_path + "plots/fit/find_Compton_peak.png")
sigmaCompton = 50

# Study of the stability of the fit changing fit extremes, in order to choose the domain of the fit
step = 30
max_step = 20
# stability_study_extreme(H, peakCompton, sigmaCompton, step, max_step, "Energy [channels]", "Counts", file_path + "plots/fit/")

min_fit = peakCompton - 5 * step
max_fit = peakCompton + 5 * step

# Study of the stability of the fit changing hist rebin, in order to choose the domain of the fit
rebin_index = [1, 2, 3, 4, 5, 6, 7, 8, 16, 32]
# stability_study_rebin(H, peakCompton, sigmaCompton, rebin_index, min_fit, max_fit, "Energy [channels]", "Counts", file_path + "plots/fit/")

# It si necessary to do the three following steps in order to have the same number of hit 
# under the peak for a rebbined histogram
hist_integral = H.Integral()
H.Rebin(3)
H.Scale(1 / 3)

fit_result, f_background, f_true = fit_peaks(H, peakCompton, sigmaCompton, min_fit, max_fit, "Energy [channels]", "Counts", file_path)

# Final fit
ll.plot_results(H, hist_integral, fit_result, f_background, f_true, min_fit, max_fit, file_path + "plots/fit/", 
                "fit_results.png", "Energy [channels]", "Counts")

