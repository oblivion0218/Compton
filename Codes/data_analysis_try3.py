import ROOT
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
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


def fit_peaks(hist, peak, sigma, min_fit, max_fit, x_axis_name, y_axis_name):
    """
    Fit a peak with a gaussian function and a background with a linear function 
    
    :param hist: ROOT histogram object.
    :param peak: Peak position.
    :param sigma: Range of the peak.
    :param min_fit: Minimum value of the fit.
    :param max_fit: Maximum value of the fit.
    :param x_axis_name: Name of the x-axis. 
    :param y_axis_name: Name of the y-axis.
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

    mpr.stampa_graph_fit(hist, f_back_e, file_path + "plots/fit/background_exp.png", "Compton peak", x_axis_name, y_axis_name, "",
                         500, 2000, 2, coo0, ["f1", "f2"])
    
    #-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    # PARTIAL FIT - Compton peak
    #-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    f_Compton = ROOT.TF1("Compton_peak", "gaus(0)", 0, 2000)
    f_Compton.SetParameter(0, 100) 
    f_Compton.SetParameter(1, peak)
    f_Compton.SetParameter(2, sigma)

    mpr.stampa_graph_fit(hist, f_Compton, file_path + "plots/fit/Compton_peak.png", "Compton peak", x_axis_name, y_axis_name, "",
                         peak - 2 * sigma, peak + 2 * sigma, 3, coo0, str0)
    
    #-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    # FIT - Complete model
    #-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    
    f_true = ROOT.TF1("model", "expo(0) + gaus(2)", 0, 2000) 
    f_true.SetParameter(0, f_back_e.GetParameter(0))
    f_true.SetParameter(1, f_back_e.GetParameter(1))
    f_true.SetParameter(2, f_Compton.GetParameter(0))
    f_true.SetParameter(3, f_Compton.GetParameter(1))
    f_true.SetParameter(4, f_Compton.GetParameter(2))

    fit_result = mpr.stampa_graph_fit(hist, f_true, file_path + "plots/fit/final_fit.png", "Spectrum", x_axis_name, y_axis_name, "",
                        min_fit, max_fit, 5, coo2, ["f1", "f2", "Amp", "<x>", "#sigma"])

    f_back_e.SetParameter(0, f_true.GetParameter(0))
    f_back_e.SetParameter(1, f_true.GetParameter(1))

    return fit_result, f_back_e, f_true


def plot_results(hist, peak, sigma, min_fit, max_fit, n_sigma_integral, fileNamePNG, x_axis_name, y_axis_name):
    """
    Plot the results of the fit."

    :param hist: ROOT histogram object."
    :param peak: Peak position."
    :param sigma: Range of the peak."
    :param min_fit: Minimum value of the fit."
    :param max_fit: Maximum value of the fit."
    :param n_sigma_integral: Range of integration in number of sigma."
    :param fileNamePNG: Name of the graph."
    :param x_axis_name: Name of the x-axis."
    :param y_axis_name: Name of the y-axis."
    """
    fit_result, f_background, f_true = fit_peaks(hist, peak, sigma, min_fit, max_fit, x_axis_name, y_axis_name)
    chi2 = fit_result.Chi2()
    ndf = fit_result.Ndf()
    chi2_ndf = chi2 / ndf

    E_mean = (f_true.GetParameter(3), f_true.GetParError(3))
    sigma = (f_true.GetParameter(4), f_true.GetParError(4))
    FWHM = (2.355 * sigma[0], 2.355 * sigma[1])
    ER = (FWHM[0] / E_mean[0], np.sqrt((FWHM[1] / E_mean[0]) ** 2 + (FWHM[0] * E_mean[1] / E_mean[0] ** 2) ** 2))

    min_integration = E_mean[0] - n_sigma_integral * sigma[0]
    max_integration = E_mean[0] + n_sigma_integral * sigma[0]
    integral = f_true.Integral(min_integration, max_integration)
    N_hit = (integral, np.sqrt(integral))
    N_hit_pc = (N_hit[0] / hist.Integral(), N_hit[1] / hist.Integral())

    text = rf"$\chi^{{2}}/\mathrm{{dof}} = {chi2:.3f}/{ndf} = {chi2_ndf:.3f}$"
    text += "\n"
    text += f"<E> = {E_mean[0]:.2f} ± {E_mean[1]:.2f}\n"
    text += f"ER = {ER[0]:.3f} ± {ER[1]:.3f}\n"
    text += f"N = {N_hit[0]:.2f} ± {N_hit[1]:.2f}\n"
    text += f"N = ({N_hit_pc[0] * 100:.3f} ± {N_hit_pc[1] * 100:.3f})%"

    #-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    # PLOT RESULTS - total fit and background
    #-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-       
    n_bins = hist.GetNbinsX()
    bin_edges = np.array([hist.GetBinLowEdge(i+1) for i in range(n_bins+1)])
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    bin_values = np.array([hist.GetBinContent(i+1) for i in range(n_bins)])

    x_values = np.linspace(min_integration, max_integration, 1000)
    f_true_func = np.vectorize(lambda x: f_true.Eval(x))
    f_back_func = np.vectorize(lambda x: f_background.Eval(x))

    y_true = f_true_func(x_values)
    y_back = f_back_func(x_values)
    
    # Create a figure with two subplots: one for the main plot and one for the residuals
    fig = plt.figure(figsize=(8, 10))
    gs = gridspec.GridSpec(2, 1, height_ratios=[3, 2])  # Main plot takes 3/4 of the space, residuals 1/4
    
    # Main plot
    ax_main = plt.subplot(gs[0])
    ax_main.hist(bin_centers, bins=bin_edges, weights=bin_values, edgecolor="gray", facecolor="none", histtype='step', label="Histogram")
    ax_main.plot(x_values, y_true, color="red", linewidth=2, label="Full Model")
    ax_main.plot(x_values, y_back, color="blue", linestyle="dashed", label="Background")
    ax_main.text(0.1, 0.75, text, fontsize=12, color="black", ha='left', transform=ax_main.transAxes)
    
    ax_main.set_xlabel(x_axis_name)
    ax_main.set_ylabel(y_axis_name)
    ax_main.legend(loc="upper right")
    ax_main.grid(True, which='both', linestyle='--', linewidth=0.5)
    ax_main.set_xticks(range(0, 2200, 200))
    
    # Residual plot
    ax_residual = plt.subplot(gs[1], sharex=ax_main)

    # Calculate residuals for points inside the integration domain
    residuals = []
    residual_centers = []
    for i, center in enumerate(bin_centers):
        if min_integration <= center <= max_integration:
            model_value = f_true.Eval(center)
            residual = (bin_values[i] - model_value) / model_value
            residuals.append(residual)
            residual_centers.append(center)  # Only include centers in the integration domain

    ax_residual.errorbar(residual_centers, residuals, fmt='x', color='black', label="Residuals")
    ax_residual.axhline(0, color="red", linewidth=2, label="Zero Line")

    ax_residual.set_xlim(min_integration - 200, max_integration + 200)
    ax_residual.set_ylabel(r"$\frac{data - model}{model}$")
    ax_residual.set_xlabel(x_axis_name)
    ax_residual.grid(True, which='both', linestyle='--', linewidth=0.5)
    ax_residual.legend(loc="upper right")

    plt.tight_layout()
    plt.savefig(file_path + "plots/fit/" + fileNamePNG)
    plt.close()


#-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# Main 
#-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
H = create_hist(file_path, "hist_sum.png")
peakCompton = ll.search_photopeak(H, 0.4, 2, file_path + "plots/fit/find_Compton_peak.png")

min_fit = []

# fit_peaks(H, peakCompton, 60, 500, 2000, "Energy [channels]", "Counts")
plot_results(H, peakCompton, 60, 900, 1300, 2, "fit_result(2sigma).png", "Energy [channels]", "Counts")

