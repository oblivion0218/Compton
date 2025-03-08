import ROOT
import numpy as np
import math
import matplotlib.pyplot as plt
from lib import MoraPyRoot as mpr
from lib import LabLibrary as ll

file_path = "/mnt/c/Users/User/Desktop/info/Compton/Measurments_trasmission/35_deg/"
coo = [0.1, 0.63, 0.45, 0.9]

# 35 --> 750 ; 50 --> 700
E_t = 750

def compute_background(E, E_t, p2, p3, p4):
    """Helper function to compute polynomial and exponential background components."""
    # Compute polynomial value and derivative at transition energy E_t
    f_pol = p2 + p3 * E_t + p4 * E_t**2
    df_pol = p3 + 2 * p4 * E_t

    if f_pol <= 0:
        return 0  # Avoid log of non-positive values

    # Compute parameters for the exponential function ensuring continuity & differentiability
    p0 = math.log(f_pol)
    p1 = df_pol / f_pol

    # Return the appropriate function value
    if E < E_t:
        return p2 + p3 * E + p4 * E**2  # Polynomial
    else:
        return math.exp(p0 + p1 * (E - E_t))  # Translated exponential


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
    peak1 = par[5] * math.exp(-((E - par[6]) ** 2) / par[7])
    peak2 = par[8] * math.exp(-((E - par[9]) ** 2) / par[10])

    return background + peak1 + peak2

# Funzione a pezzi continua e derivabile
def final_function(x, par):
    E = x[0]  # Energia (variabile indipendente)

    # Parametri del polinomio quadratico
    p2, p3, p4 = par[2], par[3], par[4]

    # Calcoliamo i parametri dell'esponenziale imponendo continuità e derivabilità
    f_pol = p2 + p3 * E_t + p4 * E_t**2
    df_pol = p3 + 2 * p4 * E_t

    # Risolviamo per p0 e p1
    p0 = np.log(f_pol)
    p1 = df_pol / f_pol

    # Applichiamo il modello corretto nelle diverse regioni
    if E < E_t:
        return p2 + p3 * E + p4 * E**2 + par[5] * np.exp(-((E - par[6]) ** 2) / par[7]) + par[8] * np.exp(-((E - par[9]) ** 2) / par[10]) + par[11] * np.exp(-((E - par[12]) ** 2) / par[13])
    else:
        return np.exp(p0 + p1 * (E - E_t)) + par[5] * np.exp(-((E - par[6]) ** 2) / par[7]) + par[8] * np.exp(-((E - par[9]) ** 2) / par[10]) + par[11] * np.exp(-((E - par[12]) ** 2) / par[13])


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


def fit_peaks(hist, peak, sigma, fileNamePNG, graph_name, x_axis_name, y_axis_name, graphic_option, pave_coordinates = None):
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
    coo0 = [0.1, 0.6, 0.4, 0.9]
    str0 = ["Amp", "<x>", "#sigma"]
    coo1 = [0.1, 0.7, 0.4, 0.9]

    n_bins = H.GetNbinsX()

    # 35 --> 0.5, 2 ; 50 --> 0.5, 2
    lead_pos =  ll.search_first_peak(hist, 0.5, 2, file_path + "plots/fit/find_lead_peak.png")
    lead_sigma = 30
    f_lead = ROOT.TF1("lead_peak", "gaus(0)", 0, 2000)
    f_lead.SetParameter(0, 0.0015)
    f_lead.SetParameter(1, lead_pos)
    f_lead.SetParameter(2, lead_sigma)

    mpr.stampa_graph_fit(hist, f_lead, file_path + "plots/fit/lead_peak.png", "Lead peak", "Energy [channels]", "Counts", "",
                         lead_pos - 2 * lead_sigma, lead_pos + 2 * lead_sigma, 3, [0.6, 0.6, 0.9, 0.9], str0)
    
    # 35 --> 0.4, 4 ; 50 --> 0.2, 4
    p511_pos = ll.search_photopeak(hist, 0.4, 4, file_path + "plots/fit/find_511_peak.png")
    p511_sigma = 50
    f_511 = ROOT.TF1("511_peak", "gaus(0)", 0, 2000)
    f_511.SetParameter(0, 0.001)
    f_511.SetParameter(1, p511_pos)
    f_511.SetParameter(2, p511_sigma)
    
    mpr.stampa_graph_fit(hist, f_511, file_path + "plots/fit/511_peak.png", "511 peak", "Energy [channels]", "Counts", "",
                         p511_pos - 2 * p511_sigma, p511_pos + 2 * p511_sigma, 3, coo0, str0)
    
    f_back_e = ROOT.TF1("background_exp", "expo(0)", 0, 2000)
    f_back_e.SetParameter(0, 1)
    f_back_e.SetParameter(1, -1)

    mpr.stampa_graph_fit(hist, f_back_e, file_path + "plots/fit/background_exp.png", "Background exp", "Energy [channels]", "Counts", "",
                        500, 2000, 2, coo1, ["p0", "p1"])
    
    f_back_p = ROOT.TF1("background_pol", "pol2(0)", 0, 2000)
    f_back_p.SetParameter(0, 1)
    f_back_p.SetParameter(1, -1)
    f_back_p.SetParameter(2, 0.001)

    mpr.stampa_graph_fit(hist, f_back_p, file_path + "plots/fit/background_pol.png", "Background pol", "Energy [channels]", "Counts", "",
                        10, 800, 3, coo1, ["a", "b", "c"])
    
    f_background = ROOT.TF1("background", background_smooth, 0, 2000, 5)
    f_background.SetParameter(0, f_back_e.GetParameter(0))
    f_background.SetParameter(1, f_back_e.GetParameter(1))
    f_background.SetParameter(2, f_back_p.GetParameter(0))
    f_background.SetParameter(3, f_back_p.GetParameter(1))
    f_background.SetParameter(4, f_back_p.GetParameter(2))

    mpr.stampa_graph_fit(hist, f_background, file_path + "plots/fit/background.png", "Background", "Energy [channels]", "Counts", "",
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

    mpr.stampa_graph_fit(hist, f_back_def, file_path + "plots/fit/background_final.png", "Background", "Energy [channels]", "Counts", "",
                        10, 2000)
    mpr.plot_TF1_MPL(f_back_def, 0, 2000, file_path + "plots/fit/background_TF1.png")

    f_Compton = ROOT.TF1("511_peak", "gaus(0)", 0, 2000)
    f_Compton.SetParameter(0, 0.001)
    f_Compton.SetParameter(1, peak)
    f_Compton.SetParameter(2, sigma)

    mpr.stampa_graph_fit(hist, f_Compton, file_path + "plots/fit/Compton_peak.png", "Compton peak", "Energy [channels]", "Counts", "",
                         peak - 2 * sigma, peak + 2 * sigma, 3, coo0, str0)
    
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

    mpr.stampa_graph_fit(hist, f_true, file_path + "plots/fit/final.png", "Spectrum", "Energy [channels]", "Counts", "",
                        10, 2000)
    mpr.plot_TF1_MPL(f_true, 0, 2000, file_path + "plots/fit/final_TF1.png")  

    for i in range(11): 
        f_back_def.SetParameter(i, f_true.GetParameter(i))

    c1 = ROOT.TCanvas("c1", "Spectrum Fit", 800, 600)

    # Draw histogram
    H.SetLineColor(ROOT.kGray + 1)  # Black color for histogram
    H.SetMarkerStyle(20)
    H.Draw("HIST SAME")  # "E" keeps error bars visible

    # Draw full spectrum function
    f_true.SetLineColor(ROOT.kRed)  # Red color for full model
    f_true.SetLineWidth(2)
    f_true.Draw("SAME")

    # Draw background function
    f_back_def.SetLineColor(ROOT.kBlue)  # Blue color for background
    f_back_def.SetLineStyle(2)  # Dashed line
    f_back_def.Draw("SAME")

    # Add legend
    legend = ROOT.TLegend(0.1, 0.65, 0.4, 0.9)
    legend.AddEntry(H, "Histogram", "lep")
    legend.AddEntry(f_back_def, "Background", "l")
    legend.AddEntry(f_true, "Full Model", "l")
    legend.Draw()

    # Save the plot
    c1.SaveAs(file_path + "plots/fit/comparison.png")

# 35 --> 0.5, 2 ; 50 --> 0.4, 2
peakCompton = ll.search_photopeak(H, 0.5, 2, file_path + "plots/fit/find_Compton_peak.png")
fit_peaks(H, peakCompton, 60, "Compton_fit.png", "Compton peak", "Energy [channels]", "Counts", "", coo)