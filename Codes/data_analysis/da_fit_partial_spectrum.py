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


# def stability_study_extreme(H, peakCompton, sigmaCompton, step, max_step, x_axis_name, y_axis_name, file_path):
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
                         file_path + "fit_stability/extreme/" + str(i))
        chi2.append(c.Chi2() / c.Ndf())

    plt.figure()
    plt.plot(n_steps, chi2, color="blue", marker="o", linestyle="dashed")
    plt.xlabel(r"$N_{step}$")
    plt.ylabel(r"$\chi^2 / ndf$")
    plt.grid()
    plt.xticks(n_steps)
    plt.savefig(file_path + "chi2_n_step.png")
    plt.close()


def stability_study_extreme(H, peakCompton, sigmaCompton, step, max_step, x_axis_name, y_axis_name, file_path):
    """
    Study the stability of the fit as a function of the fitting range extremes.
    In addition to χ²/ndf, also evaluate the fitted Compton peak position, σ, and
    the energy resolution (ER = σ/peak × 100%).
    
    :param H: ROOT histogram object.
    :param peakCompton: Nominal peak position.
    :param sigmaCompton: Nominal sigma of the peak.
    :param step: Step size for changing the fit extremes.
    :param max_step: Maximum step (number of iterations).
    :param x_axis_name: Label for x-axis.
    :param y_axis_name: Label for y-axis.
    :param file_path: Path (folder) to save the plots.
    """
    n_steps = [i for i in range(1, max_step)]
    chi2_list = []
    peak_fit_list = []
    sigma_fit_list = []
    er_list = []
    
    # Loop over fit extremes
    for i in n_steps:
        min_fit = peakCompton - i * step
        max_fit = peakCompton + i * step
        # Call your fit_peaks routine (which returns a canvas, background fit, and the full fit function)
        c, _, f_true = fit_peaks(H, peakCompton, sigmaCompton,
                                 min_fit, max_fit,
                                 x_axis_name, y_axis_name,
                                 file_path + "fit_stability/extreme/" + str(i))
        chi2_val = c.Chi2() / c.Ndf()
        fitted_peak = f_true.GetParameter(3)  # Compton peak position
        fitted_sigma = f_true.GetParameter(4)  # Compton sigma
        er = (fitted_sigma * 2.355 / fitted_peak) * 100  # Energy resolution in %
        
        chi2_list.append(chi2_val)
        peak_fit_list.append(fitted_peak)
        sigma_fit_list.append(fitted_sigma)
        er_list.append(er)
    
    # Plot all the stability studies in a single figure using subplots.
    fig, axs = plt.subplots(2, 2, figsize=(18, 10))
    
    # Subplot 1: χ²/ndf stability
    axs[0, 0].plot(n_steps, chi2_list, marker='o', linestyle='--', color='blue')
    axs[0, 0].set_title(r'$\chi^2/ndf$ Stability')
    axs[0, 0].set_xlabel("N_step")
    axs[0, 0].set_ylabel(r"$\chi^2/ndf$")
    axs[0, 0].grid(True)
    axs[0, 0].set_xticks(n_steps)
    
    # Subplot 2: Fitted Peak Position
    axs[0, 1].plot(n_steps, peak_fit_list, marker='o', linestyle='--', color='red')
    axs[0, 1].set_title("Fitted Compton Peak Position")
    axs[0, 1].set_xlabel("N_step")
    axs[0, 1].set_ylabel("Peak Position")
    axs[0, 1].grid(True)
    axs[0, 1].set_xticks(n_steps)
    
    # Subplot 3: Fitted Sigma
    axs[1, 0].plot(n_steps, sigma_fit_list, marker='o', linestyle='--', color='green')
    axs[1, 0].set_title("Fitted σ Stability")
    axs[1, 0].set_xlabel("N_step")
    axs[1, 0].set_ylabel("σ")
    axs[1, 0].grid(True)
    axs[1, 0].set_xticks(n_steps)
    
    # Subplot 4: Energy Resolution (ER = σ/peak × 100)
    axs[1, 1].plot(n_steps, er_list, marker='o', linestyle='--', color='magenta')
    axs[1, 1].set_title("Energy Resolution Stability")
    axs[1, 1].set_xlabel("N_step")
    axs[1, 1].set_ylabel("ER (%)")
    axs[1, 1].grid(True)
    axs[1, 1].set_xticks(n_steps)
    
    fig.tight_layout()
    plot_file = file_path + "stability_summary_extreme.png"
    plt.savefig(plot_file)
    plt.close()
    print(f"Stability summary plot saved to {plot_file}")


def stability_study_rebin(H, peakCompton, sigmaCompton, rebin_max, min_fit, max_fit, x_axis_name, y_axis_name, file_path):
    """
    Study the stability of the fit as a function of the histogram rebinning factor.
    In addition to χ²/ndf, evaluate the fitted Compton peak position, σ, and the energy 
    resolution (ER = σ/peak × 100%).
    
    :param H: ROOT histogram object.
    :param peakCompton: Nominal peak position.
    :param sigmaCompton: Nominal sigma of the peak.
    :param rebin_max: Maximum rebin factor to loop over (from 1 to rebin_max).
    :param min_fit: Minimum value of the fit.
    :param max_fit: Maximum value of the fit.
    :param x_axis_name: Label for x-axis.
    :param y_axis_name: Label for y-axis.
    :param file_path: Path (folder) to save the plots.
    """

    # Create lists to store stability parameters
    rebin_factors = [i for i in range(1, rebin_max + 1)]
    chi2_list = []
    peak_fit_list = []
    sigma_fit_list = []
    er_list = []
    
    # Loop over rebin factors
    for i in rebin_factors:
        # Clone and rebin the histogram to avoid modifying the original
        H_rebin = H.Clone()
        H_rebin.Rebin(i)
        
        # Perform the fit on the rebinned histogram
        c, _, f_true = fit_peaks(H_rebin, peakCompton, sigmaCompton,
                                 min_fit, max_fit,
                                 x_axis_name, y_axis_name,
                                 file_path + "fit_stability/rebin/" + str(i))
        # Extract χ²/ndf
        chi2_val = c.Chi2() / c.Ndf()
        # Extract fitted parameters (using indices as in your fit_peaks)
        fitted_peak = f_true.GetParameter(3)  # Compton peak position
        fitted_sigma = f_true.GetParameter(4)   # Compton sigma
        er = (fitted_sigma * 2.355 / fitted_peak) * 100   # Energy resolution in %
        
        chi2_list.append(chi2_val)
        peak_fit_list.append(fitted_peak)
        sigma_fit_list.append(fitted_sigma)
        er_list.append(er)
    
    # Plot all the stability studies in a single figure using subplots
    fig, axs = plt.subplots(2, 2, figsize=(18, 10))
    
    # Subplot 1: χ²/ndf stability vs. rebin factor
    axs[0, 0].plot(rebin_factors, chi2_list, marker='o', linestyle='--', color='blue')
    axs[0, 0].set_title(r'$\chi^2/ndf$ Stability vs Rebin Factor')
    axs[0, 0].set_xlabel("Rebin Factor")
    axs[0, 0].set_ylabel(r"$\chi^2/ndf$")
    axs[0, 0].grid(True)
    axs[0, 0].set_xticks(rebin_factors)
    
    # Subplot 2: Fitted Peak Position vs. rebin factor
    axs[0, 1].plot(rebin_factors, peak_fit_list, marker='o', linestyle='--', color='red')
    axs[0, 1].set_title("Fitted Compton Peak Position vs Rebin Factor")
    axs[0, 1].set_xlabel("Rebin Factor")
    axs[0, 1].set_ylabel("Peak Position")
    axs[0, 1].grid(True)
    axs[0, 1].set_xticks(rebin_factors)
    
    # Subplot 3: Fitted σ vs. rebin factor
    axs[1, 0].plot(rebin_factors, sigma_fit_list, marker='o', linestyle='--', color='green')
    axs[1, 0].set_title("Fitted σ Stability vs Rebin Factor")
    axs[1, 0].set_xlabel("Rebin Factor")
    axs[1, 0].set_ylabel("σ")
    axs[1, 0].grid(True)
    axs[1, 0].set_xticks(rebin_factors)
    
    # Subplot 4: Energy Resolution vs. rebin factor
    axs[1, 1].plot(rebin_factors, er_list, marker='o', linestyle='--', color='magenta')
    axs[1, 1].set_title("Energy Resolution vs Rebin Factor")
    axs[1, 1].set_xlabel("Rebin Factor")
    axs[1, 1].set_ylabel("ER (%)")
    axs[1, 1].grid(True)
    axs[1, 1].set_xticks(rebin_factors)
    
    fig.tight_layout()
    plot_file = file_path + "stability_summary_rebin.png"
    plt.savefig(plot_file)
    plt.close()
    print(f"Rebin stability summary plot saved to {plot_file}")


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

n_steps = 5
min_fit = peakCompton - n_steps * step
max_fit = peakCompton + n_steps * step

# Study of the stability of the fit changing hist rebin, in order to choose the domain of the fit
rebin_max = 32
# stability_study_rebin(H, peakCompton, sigmaCompton, rebin_max, min_fit, max_fit, "Energy [channels]", "Counts", file_path + "plots/fit/")

hist_integral = H.Integral()
rebin_param = 5
H.Rebin(rebin_param)

fit_result, f_background, f_true = fit_peaks(H, peakCompton, sigmaCompton, min_fit, max_fit, "Energy [channels]", "Counts", 
                                             file_path + "plots/fit/")

# Final fit
ll.plot_results(H, hist_integral, fit_result, f_background, f_true, rebin_param, min_fit, max_fit, file_path + "plots/fit/", 
                "fit_results.png", "Energy [channels]", "Counts")

