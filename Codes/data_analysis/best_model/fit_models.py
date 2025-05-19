import ROOT
from lib import MoraPyRoot as mpr
from lib import LabLibrary as ll

# FIT MODELS
# exp + gaus 
# p1 + gaus
# p2 + gaus
# p4 + gaus

def fit_exp_back(hist, amp, peak, sigma, min_fit, max_fit, x_axis_name, y_axis_name, file_path):
    f_background = ROOT.TF1("background", "expo(0)", 0, 2000)
    f_background.SetParameter(0, 1)
    f_background.SetParameter(1, 1)

    f_peak = ROOT.TF1("peak", "gaus(0)", 0, 2000)
    f_peak.SetParameter(0, amp)
    f_peak.SetParameter(1, peak)
    f_peak.SetParameter(2, sigma)

    f_true = ROOT.TF1("model", "expo(0) + gaus(2)", 0, 2000) 
    f_true.SetParameter(0, f_background.GetParameter(0))
    f_true.SetParameter(1, f_background.GetParameter(1))
    f_true.SetParameter(2, f_peak.GetParameter(0))
    f_true.SetParameter(3, f_peak.GetParameter(1))
    f_true.SetParameter(4, f_peak.GetParameter(2))

    fit_result = mpr.stampa_graph_fit(hist, f_true, file_path + "fit_exp.png", "Spectrum", 
                                      x_axis_name, y_axis_name, "", min_fit, max_fit)

    f_background.SetParameter(0, f_true.GetParameter(0))
    f_background.SetParameter(1, f_true.GetParameter(1))

    return fit_result, f_background, f_true

def fit_p1_back(hist, amp, peak, sigma, min_fit, max_fit, x_axis_name, y_axis_name, file_path):
    f_background = ROOT.TF1("background", "pol1(0)", 0, 2000)
    f_background.SetParameter(0, 1)
    f_background.SetParameter(1, 1)

    f_peak = ROOT.TF1("peak", "gaus(0)", 0, 2000)
    f_peak.SetParameter(0, amp)
    f_peak.SetParameter(1, peak)
    f_peak.SetParameter(2, sigma)

    f_true = ROOT.TF1("model", "pol1(0) + gaus(2)", 0, 2000) 
    f_true.SetParameter(0, f_background.GetParameter(0))
    f_true.SetParameter(1, f_background.GetParameter(1))
    f_true.SetParameter(2, f_peak.GetParameter(0))
    f_true.SetParameter(3, f_peak.GetParameter(1))
    f_true.SetParameter(4, f_peak.GetParameter(2))

    fit_result = mpr.stampa_graph_fit(hist, f_true, file_path + "fit_p1.png", "Spectrum", 
                                      x_axis_name, y_axis_name, "", min_fit, max_fit)

    f_background.SetParameter(0, f_true.GetParameter(0))
    f_background.SetParameter(1, f_true.GetParameter(1))

    return fit_result, f_background, f_true

def fit_p2_back(hist, amp, peak, sigma, min_fit, max_fit, x_axis_name, y_axis_name, file_path):
    f_background = ROOT.TF1("background", "pol2(0)", 0, 2000)
    f_background.SetParameter(0, 1)
    f_background.SetParameter(1, 1)
    f_background.SetParameter(2, 1)

    f_peak = ROOT.TF1("peak", "gaus(0)", 0, 2000)
    f_peak.SetParameter(0, amp)
    f_peak.SetParameter(1, peak)
    f_peak.SetParameter(2, sigma)

    f_true = ROOT.TF1("model", "pol2(0) + gaus(3)", 0, 2000) 
    f_true.SetParameter(0, f_background.GetParameter(0))
    f_true.SetParameter(1, f_background.GetParameter(1))
    f_true.SetParameter(2, f_background.GetParameter(2))
    f_true.SetParameter(3, f_peak.GetParameter(0))
    f_true.SetParameter(4, f_peak.GetParameter(1))
    f_true.SetParameter(5, f_peak.GetParameter(2))

    fit_result = mpr.stampa_graph_fit(hist, f_true, file_path + "fit_p2.png", "Spectrum", 
                                      x_axis_name, y_axis_name, "", min_fit, max_fit)

    f_background.SetParameter(0, f_true.GetParameter(0))
    f_background.SetParameter(1, f_true.GetParameter(1))
    f_background.SetParameter(2, f_true.GetParameter(2))

    return fit_result, f_background, f_true

def fit_p4_back(hist, amp, peak, sigma, min_fit, max_fit, x_axis_name, y_axis_name, file_path):
    f_background = ROOT.TF1("background", "pol4(0)", 0, 2000)
    f_background.SetParameter(0, 1)
    f_background.SetParameter(1, 1)
    f_background.SetParameter(2, 1)
    f_background.SetParameter(3, 1)

    f_peak = ROOT.TF1("peak", "gaus(0)", 0, 2000)
    f_peak.SetParameter(0, amp)
    f_peak.SetParameter(1, peak)
    f_peak.SetParameter(2, sigma)

    f_true = ROOT.TF1("model", "pol4(0) + gaus(4)", 0, 2000) 
    f_true.SetParameter(0, f_background.GetParameter(0))
    f_true.SetParameter(1, f_background.GetParameter(1))
    f_true.SetParameter(2, f_background.GetParameter(2))
    f_true.SetParameter(3, f_background.GetParameter(3))
    f_true.SetParameter(4, f_peak.GetParameter(0))
    f_true.SetParameter(5, f_peak.GetParameter(1))
    f_true.SetParameter(6, f_peak.GetParameter(2))

    fit_result = mpr.stampa_graph_fit(hist, f_true, file_path + "fit_p4.png", "Spectrum", 
                                      x_axis_name, y_axis_name, "", min_fit, max_fit)

    f_background.SetParameter(0, f_true.GetParameter(0))
    f_background.SetParameter(1, f_true.GetParameter(1))
    f_background.SetParameter(2, f_true.GetParameter(2))
    f_background.SetParameter(3, f_true.GetParameter(3))

    return fit_result, f_background, f_true


