import ROOT
import numpy as np
import matplotlib.pyplot as plt
from lib import MoraPyRoot as mpr
from lib import LabLibrary as ll

#-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# Parameters to set for begin
#-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
file_path = "/mnt/c/Users/User/Desktop/info/Compton/Measurments_trasmission/50_deg/"

# E_t = 750 # 35°
E_t = 700 # 50°

# [lead peak , 511 peak, Compton peak]
# find_peak_param = [(0.5, 2), (0.4, 4), (0.5, 2)] # 35°
find_peak_param = [(0.5, 2), (0.2, 4), (0.4, 2)] # 50°

def compute_background(E, E_t, p2, p3, p4):
    """Helper function to compute polynomial and exponential background components."""
    # Compute polynomial value and derivative at transition energy E_t
    f_pol = p2 + p3 * E_t + p4 * E_t**2
    df_pol = p3 + 2 * p4 * E_t

    if f_pol <= 0:
        return 0  # Avoid log of non-positive values

    # Compute parameters for the exponential function ensuring continuity & differentiability
    p0 = np.log(f_pol)
    p1 = df_pol / f_pol

    # Return the appropriate function value
    if E < E_t:
        return p2 + p3 * E + p4 * E**2  # Polynomial
    else:
        return np.exp(p0 + p1 * (E - E_t))  # Translated exponential


def background_smooth(x, par):
    """Background function without peaks."""
    E = x[0]  # Independent variable (energy)

    # Quadratic polynomial parameters
    p2, p3, p4 = par[2], par[3], par[4]

    return compute_background(E, E_t, p2, p3, p4)


def background_smooth_peak(x, par):
    """Background function with additional Gaussian peaks."""
    E = x[0]  # Independent variable (energy)

    # Quadratic polynomial parameters
    p2, p3, p4 = par[2], par[3], par[4]

    # Compute base background function
    background = compute_background(E, E_t, p2, p3, p4)

    # Add Gaussian peaks
    peak1 = par[5] * np.exp(-((E - par[6]) ** 2) / par[7])
    peak2 = par[8] * np.exp(-((E - par[9]) ** 2) / par[10])

    return background + peak1 + peak2


def final_function(x, par):
    """Background function with additional Gaussian peaks."""
    E = x[0]  # Independent variable (energy)

    # Quadratic polynomial parameters
    p2, p3, p4 = par[2], par[3], par[4]

    # Compute base background function
    background = compute_background(E, E_t, p2, p3, p4)

    # Add Gaussian peaks
    peak1 = par[5] * np.exp(-((E - par[6]) ** 2) / par[7])
    peak2 = par[8] * np.exp(-((E - par[9]) ** 2) / par[10])
    peak3 = par[11] * np.exp(-((E - par[12]) ** 2) / par[13])

    return background + peak1 + peak2 + peak3


def create_hist(file_path, fileNamePNG):
    hs = ll.hist_vector(file_path + "data/")

    for h in hs:
        mpr.plot_hist_MPL(h, file_path + "plots/hist/" + h.GetName() + ".png")
    
    H = ll.spectum_sum(hs)
    mpr.plot_hist_MPL(H, file_path + "plots/hist/" + fileNamePNG)

    return H


H = create_hist(file_path, "hist_sum.png")
integral_H = H.Integral()
H = ll.normalize_histogram(H)


def fit_peaks(hist, peak, sigma, n_sigma_integral, fileNamePNG, x_axis_name, y_axis_name):
    """
    Fit a peak with a gaussian function and a background with a linear function 
    
    :param hist: ROOT histogram object.
    :param peak: Peak position.
    :param sigma: Range of the peak.
    :param n_sigma_integral: Range of integration in number of sigma.
    :param fileNamePNG: Name of the graph.
    :param x_axis_name: Name of the x-axis. 
    :param y_axis_name: Name of the y-axis.
    """
    coo0 = [0.1, 0.6, 0.4, 0.9]
    str0 = ["Amp", "<x>", "#sigma"]
    coo1 = [0.1, 0.7, 0.4, 0.9]

    #-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    # PARTIAL FIT - Background peaks
    #-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    lead_pos =  ll.search_first_peak(hist, find_peak_param[0][0], find_peak_param[0][1], file_path + "plots/fit/find_lead_peak.png")
    lead_sigma = 30
    f_lead = ROOT.TF1("lead_peak", "gaus(0)", 0, 2000)
    f_lead.SetParameter(0, 0.0015)
    f_lead.SetParameter(1, lead_pos)
    f_lead.SetParameter(2, lead_sigma)

    mpr.stampa_graph_fit(hist, f_lead, file_path + "plots/fit/lead_peak.png", "Lead peak", x_axis_name, y_axis_name, "",
                         lead_pos - 2 * lead_sigma, lead_pos + 2 * lead_sigma, 3, [0.6, 0.6, 0.9, 0.9], str0)
    
    p511_pos = ll.search_photopeak(hist, find_peak_param[1][0], find_peak_param[1][1], file_path + "plots/fit/find_511_peak.png")
    p511_sigma = 50
    f_511 = ROOT.TF1("511_peak", "gaus(0)", 0, 2000)
    f_511.SetParameter(0, 0.001)
    f_511.SetParameter(1, p511_pos)
    f_511.SetParameter(2, p511_sigma)
    
    mpr.stampa_graph_fit(hist, f_511, file_path + "plots/fit/511_peak.png", "511 peak", x_axis_name, y_axis_name, "",
                         p511_pos - 2 * p511_sigma, p511_pos + 2 * p511_sigma, 3, coo0, str0)
    
    #-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    # PARTIAL FIT - Background
    #-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    f_back_e = ROOT.TF1("background_exp", "expo(0)", 0, 2000)
    f_back_e.SetParameter(0, 1)
    f_back_e.SetParameter(1, -1)

    mpr.stampa_graph_fit(hist, f_back_e, file_path + "plots/fit/background_exp.png", "Background exp", x_axis_name, y_axis_name, "",
                        500, 2000, 2, coo1, ["p0", "p1"])
    
    f_back_p = ROOT.TF1("background_pol", "pol2(0)", 0, 2000)
    f_back_p.SetParameter(0, 1)
    f_back_p.SetParameter(1, -1)
    f_back_p.SetParameter(2, 0.001)

    mpr.stampa_graph_fit(hist, f_back_p, file_path + "plots/fit/background_pol.png", "Background pol", x_axis_name, y_axis_name, "",
                        10, 800, 3, coo1, ["a", "b", "c"])
    
    f_background = ROOT.TF1("background", background_smooth, 0, 2000, 5)
    f_background.SetParameter(0, f_back_e.GetParameter(0))
    f_background.SetParameter(1, f_back_e.GetParameter(1))
    f_background.SetParameter(2, f_back_p.GetParameter(0))
    f_background.SetParameter(3, f_back_p.GetParameter(1))
    f_background.SetParameter(4, f_back_p.GetParameter(2))

    mpr.stampa_graph_fit(hist, f_background, file_path + "plots/fit/background.png", "Background", x_axis_name, y_axis_name, "",
                        10, 2000)
 
    f_back_def = ROOT.TF1("final background", background_smooth_peak, 0, 2000, 11)
    f_back_def.SetParameter(0, f_background.GetParameter(0))
    f_back_def.SetParameter(1, f_background.GetParameter(1))
    f_back_def.SetParameter(2, f_background.GetParameter(2))
    f_back_def.SetParameter(3, f_background.GetParameter(3))
    f_back_def.SetParameter(4, f_background.GetParameter(4))
    f_back_def.SetParameter(5, f_lead.GetParameter(0))
    f_back_def.SetParameter(6, f_lead.GetParameter(1))
    f_back_def.SetParameter(7, f_lead.GetParameter(2))
    f_back_def.SetParameter(8, f_511.GetParameter(0))
    f_back_def.SetParameter(9, f_511.GetParameter(1))
    f_back_def.SetParameter(10, f_511.GetParameter(2))

    mpr.stampa_graph_fit(hist, f_back_def, file_path + "plots/fit/background_final.png", "Background", x_axis_name, y_axis_name, "",
                        10, 2000)

    #-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    # PARTIAL FIT - Compton peak
    #-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    f_Compton = ROOT.TF1("Compton_peak", "gaus(0)", 0, 2000)
    f_Compton.SetParameter(0, 0.001)
    f_Compton.SetParameter(1, peak)
    f_Compton.SetParameter(2, sigma)

    mpr.stampa_graph_fit(hist, f_Compton, file_path + "plots/fit/Compton_peak.png", "Compton peak", x_axis_name, y_axis_name, "",
                         peak - 2 * sigma, peak + 2 * sigma, 3, coo0, str0)
    
    #-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    # FIT - Complete model
    #-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    f_true = ROOT.TF1("model", final_function, 0, 2000, 14)
    f_true.SetParameter(0, f_back_def.GetParameter(0))
    f_true.SetParameter(1, f_back_def.GetParameter(1))
    f_true.SetParameter(2, f_back_def.GetParameter(2))
    f_true.SetParameter(3, f_back_def.GetParameter(3))
    f_true.SetParameter(4, f_back_def.GetParameter(4))
    f_true.SetParameter(5, f_back_def.GetParameter(5))
    f_true.SetParameter(6, f_back_def.GetParameter(6))
    f_true.SetParameter(7, f_back_def.GetParameter(7))
    f_true.SetParameter(8, f_back_def.GetParameter(8))
    f_true.SetParameter(9, f_back_def.GetParameter(9))
    f_true.SetParameter(10, f_back_def.GetParameter(10))
    f_true.SetParameter(11, f_Compton.GetParameter(0))
    f_true.SetParameter(12, f_Compton.GetParameter(1))
    f_true.SetParameter(13, f_Compton.GetParameter(2))

    mpr.stampa_graph_fit(hist, f_true, file_path + "plots/fit/final_fit.png", "Spectrum", x_axis_name, y_axis_name, "",
                        10, 2000)
    mpr.plot_TF1_MPL(f_true, 0, 2000, file_path + "plots/fit/final_TF1.png")  

    for i in range(11): 
        f_back_def.SetParameter(i, f_true.GetParameter(i))
    
    mpr.plot_TF1_MPL(f_back_def, 0, 2000, file_path + "plots/fit/background_TF1.png")

    #-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    # Relevant parameters
    #-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    E_mean = (f_true.GetParameter(12), f_true.GetParError(12))
    FWHM = (2.355 * f_true.GetParameter(13), 2.355 * f_true.GetParError(13))
    ER = (FWHM[0] / E_mean[0], np.sqrt((FWHM[1] / E_mean[0]) ** 2 + (FWHM[0] * E_mean[1] / E_mean[0] ** 2) ** 2))
    min = E_mean[0] - n_sigma_integral * f_true.GetParameter(13)
    max = E_mean[0] + n_sigma_integral * f_true.GetParameter(13)
    N_pc = f_true.Integral(min, max) - f_back_def.Integral(min, max)
    N = integral_H * N_pc
    N_hit = (N, np.sqrt(N))
    N_hit_pc = (N_pc, N_hit[1] / integral_H)

    text = f"<E> = {E_mean[0]:.0f} ± {E_mean[1]:.0f}\n"
    text += f"ER = {ER[0]:.3f} ± {ER[1]:.3f}\n"
    text += f"N = {N_hit[0]:.0f} ± {N_hit[1]:.0f}\n"
    text += f"N = ({N_hit_pc[0] * 100:.3f} ± {N_hit_pc[1] * 100:.3f})%"

    #-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    # Print result
    #-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    bin_edges = np.array([H.GetBinLowEdge(i+1) for i in range(H.GetNbinsX()+1)])
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    bin_values = np.array([H.GetBinContent(i+1) for i in range(H.GetNbinsX())])

    # Define x-range for functions
    x_values = np.linspace(bin_edges[0], bin_edges[-1], 1000)

    # Assuming f_true and f_back_def are TF1 (convert them to callable functions)
    f_true_func = np.vectorize(lambda x: f_true.Eval(x))
    f_back_func = np.vectorize(lambda x: f_back_def.Eval(x))

    # Compute function values
    y_true = f_true_func(x_values)
    y_back = f_back_func(x_values)

    # Create figure
    plt.figure(figsize=(8, 6))

    # Plot histogram
    plt.hist(bin_centers, bins=bin_edges, weights=bin_values, edgecolor="gray", facecolor="none", histtype='step', label="Histogram")

    # Plot functions
    plt.plot(x_values, y_true, color="red", linewidth=2, label="Full Model")
    plt.plot(x_values, y_back, color="blue", linestyle="dashed", label="Background")

    # Add the text to the plot
    plt.text(0.2, 0.8, text, fontsize=12, color="black", ha='left', transform=plt.gca().transAxes)

    # Labels and title
    plt.xlabel(x_axis_name)
    plt.ylabel(y_axis_name)

    # Legend
    plt.legend(loc="upper right")

    # Save plot
    plt.savefig(file_path + "plots/fit/" + fileNamePNG)


peakCompton = ll.search_photopeak(H, find_peak_param[2][0], find_peak_param[2][1], file_path + "plots/fit/find_Compton_peak.png")
fit_peaks(H, peakCompton, 60, 2, "fit_result(2sigma).png", "Energy [channels]", "Counts")