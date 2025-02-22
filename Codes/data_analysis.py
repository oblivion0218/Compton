import ROOT
import numpy as np
import math
from lib import MoraPyRoot as mpr
from lib import LabLibrary as ll

file_path = "/mnt/c/Users/User/Desktop/info/Compton/Measurments_trasmission/50_deg/"

hs = ll.hist_vector(file_path + "data/")

for h in hs:
    mpr.plot_hist_MPL(h, file_path + "plots/hist/" + h.GetName() + ".png")
H = ll.spectum_sum(hs)
mpr.plot_hist_MPL(H, file_path + "plots/hist/" + "hist_sum.png")

def fit_peaks(hist, peak, extreme, fileNamePNG, graph_name, x_axis_name, y_axis_name, graphic_option, pave_coordinates = None):
    coo2 = [0.1, 0.5, 0.4, 0.9]
    str2 = ["f1", "f2", "Amp", "<x>", "#sigma"]

    n_bins = H.GetNbinsX()

    f_picco = ROOT.TF1("picco", "gaus(0)", 0, 2500)
    # L'inizializzazione dei parametri varia al variare del numero di bin
    f_picco.SetParameters(0, 40000 / n_bins)
    f_picco.SetParameters(1, peak)
    f_picco.SetParameters(2, 0.001 * n_bins)
    mpr.stampa_graph_fit(H, f_picco, file_path + "plots/peak_" + fileNamePNG, "", "", "Counts", "", peak - extreme, peak + extreme)

    f_fondo = ROOT.TF1("fondo", "[0] + [1]/x", 0, 2500)
    f_fondo.SetParameters(0, 1)
    f_fondo.SetParameters(1, 1)
    mpr.stampa_graph_fit(H, f_fondo, file_path + "plots/background_" + fileNamePNG, "", "", "Counts", "", peak - extreme * 3, peak + extreme * 3)

    # Background + peak
    f_true = ROOT.TF1("modello", "[0] + [1]/x + gaus(2)", 0, 2500)
    f_true.SetParameter(0, f_fondo.GetParameter(0))
    f_true.SetParameter(1, f_fondo.GetParameter(1))
    f_true.SetParameter(2, f_picco.GetParameter(0))
    f_true.SetParameter(3, f_picco.GetParameter(1))
    f_true.SetParameter(4, f_picco.GetParameter(2))

    canvas = ROOT.TCanvas()
    
    hist.Draw(graphic_option)
    hist.SetTitle(graph_name)
    hist.GetXaxis().SetTitle(x_axis_name)
    hist.GetYaxis().SetTitle(y_axis_name)

    min_val = peak - extreme
    max_val = peak + extreme

    fit_result = hist.Fit(f_true, "S", "", min_val, max_val) 

    if pave_coordinates != None: 
        text_box = ROOT.TPaveText(pave_coordinates[0], pave_coordinates[1], pave_coordinates[2], pave_coordinates[3], "NDC")
        text_box.SetFillColor(0)
        text_box.SetTextAlign(12)

        E_mean = f_true.GetParameter(3)
        E_mean_error = f_true.GetParError(3)
        E_mean_SN = mpr.exponential(E_mean_error)

        if abs(E_mean_SN.exp) < 3:
            text_box.AddText(f"<E> = {E_mean:.3f} #pm {E_mean_error:.3f}")
        else:   
            text_box.AddText(f"E_mean = ({E_mean * 10 ** -E_mean_SN.exp:.3f} #pm {E_mean_SN.n:.3f}) * 10^{{{E_mean_SN.exp}}}")

        FWHM = 2.355 * f_true.GetParameter(4)

        ER = FWHM / E_mean
        ER_error = np.sqrt((2.355 * f_true.GetParError(4) / E_mean)**2 + (FWHM * f_true.GetParError(3) / E_mean**2)**2)
        ER_SN = mpr.exponential(ER_error)

        if abs(ER_SN.exp) < 3:
            text_box.AddText(f"ER = {ER:.3f} #pm {ER_error:.3f}")
        else:    
            text_box.AddText(f"ER = ({ER * 10 ** -ER_SN.exp:.3f} #pm {ER_SN.n:.3f}) * 10^{{{ER_SN.exp}}}")

        N_hit = f_true.Integral(peak - extreme, peak + extreme) - f_fondo.Integral(peak - extreme, peak + extreme)
        N_hit_error = np.sqrt(N_hit)
        N_hit_SN = mpr.exponential(N_hit_error)

        if abs(N_hit_SN.exp) < 3:
            text_box.AddText(f"N = {N_hit:.0f} #pm {N_hit_error:.3f}")
        else:
            text_box.AddText(f"N = ({N_hit * 10 ** -N_hit_SN.exp:.3f} #pm {N_hit_SN.n:.3f}) * 10^{{{N_hit_SN.exp}}}")

        text_box.Draw()

        canvas.Print(file_path + "plots/" + fileNamePNG, "png")
        del text_box
    
    else:            
        canvas.Print(file_path + "plots/" + fileNamePNG, "png")

peakCompton = ll.search_photopeak(H, 0.4, 4)
fit_peaks(H, peakCompton, 150, "Compton_fit.png", "Compton peak", "Energy [keV]", "Counts", "", [0.1, 0.68, 0.4, 0.9])

peak511 = ll.search_photopeak(H, 0.2, 4)
fit_peaks(H, peak511, 100, "511_fit.png", "511 peak", "Energy [keV]", "Counts", "", [0.1, 0.68, 0.4, 0.9])
