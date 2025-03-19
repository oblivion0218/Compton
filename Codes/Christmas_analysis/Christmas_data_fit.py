import ROOT
import pandas as pd
import numpy as np
from lib import MoraPyRoot as mpr
from lib import LabLibrary as ll

#Andrea
file_path = "/mnt/c/Users/User/Desktop/info/Compton/Christmas_measurments/"

data = pd.read_csv(file_path + "results/Misure_gennaio_3.txt", delim_whitespace=True)

RE = data["RE"]
RE_error = data["sigma_RE"]
E_mean = data["<E>"]
E_mean_error = data["sigma_<E>"]

times = []
s = 0
for i in range(len(RE)):
    minutes = 20

    s += minutes
    times.append(s)

T_RE = ROOT.TGraphErrors(len(times), np.array(times, dtype='d'), np.array(RE, dtype='d'), np.zeros(len(times), dtype='d'), np.array(RE_error, dtype='d'))
T_E = ROOT.TGraphErrors(len(times), np.array(times, dtype='d'), np.array(E_mean, dtype='d'), np.zeros(len(times), dtype='d'), np.array(E_mean_error, dtype='d'))

f_trueRE = ROOT.TF1("f_trueRE", "[0] + [1] * x", 2500, 8000)
f_trueRE.SetParameter(0, 0.06875)
f_trueRE.FixParameter(1, 0)

mpr.fit(T_RE, f_trueRE, 1, "S", 1000, 2500, 8000)

f_trueE = ROOT.TF1("f_trueE", "[0] + [1] * x", 2500, 8000)
f_trueE.SetParameter(0, 1430)
f_trueE.FixParameter(1, 0)

f_trueE.SetParameter(0, 1430)

mpr.fit(T_E, f_trueE, 1, "S", 1000, 2500, 8000)

def plot_fit_results(T_RE, T_E, f_RE, f_E):
    """
    Plots the TGraphErrors and their respective fits with legends.
    Prints the fit parameters to the console.

    Parameters:
        T_RE (ROOT.TGraphErrors): Graph of RE vs time.
        T_E (ROOT.TGraphErrors): Graph of E_mean vs time.
        f_RE (ROOT.TF1): Fit function for RE.
        f_E (ROOT.TF1): Fit function for E_mean.
    """
    # Create a canvas with two pads
    canvas = ROOT.TCanvas("canvas", "Fit Results", 1200, 800)
    canvas.Divide(1, 2)

    # Plot T_RE on the first pad
    canvas.cd(1)
    ROOT.gPad.SetGrid()
    T_RE.SetTitle("Energetic Resolution vs Time;Time (s);ER")
    T_RE.SetMarkerStyle(20)
    T_RE.SetMarkerColor(ROOT.kBlue)
    T_RE.Draw("AP")
    f_RE.SetLineColor(ROOT.kRed)
    f_RE.SetLineWidth(2)
    f_RE.Draw("SAME")

    # Add a legend for T_RE
    legend_RE = ROOT.TLegend(0.6, 0.65, 0.9, 0.9)
    legend_RE.SetHeader("Fit Parameters for ER", "C")
    legend_RE.AddEntry(T_RE, "Data Points", "p")
    legend_RE.AddEntry(f_RE, f"Fit: q = {f_RE.GetParameter(0):.6f} #pm {f_RE.GetParError(0):.6f}", "l")
    legend_RE.Draw()

    # Plot T_E on the second pad
    canvas.cd(2)
    ROOT.gPad.SetGrid()
    T_E.SetTitle("Mean Energy vs Time;Time (s);<E> (keV)")
    T_E.SetMarkerStyle(20)
    T_E.SetMarkerColor(ROOT.kGreen)
    T_E.Draw("AP")
    f_E.SetLineColor(ROOT.kRed)
    f_E.SetLineWidth(2)
    f_E.Draw("SAME")

    # Add a legend for T_E
    legend_E = ROOT.TLegend(0.6, 0.65, 0.9, 0.9)
    legend_E.SetHeader("Fit Parameters for <E>", "C")
    legend_E.AddEntry(T_E, "Data Points", "p")
    legend_E.AddEntry(f_E, f"Fit: q = {f_E.GetParameter(0):.4f} #pm {f_E.GetParError(0):.4f}", "l")
    legend_E.Draw()

    # Update and show the canvas
    canvas.Update()
    canvas.SaveAs(file_path + "plots/fit_results.png")  # Optionally save the canvas as an image
    canvas.Show()

# Call the function
plot_fit_results(T_RE, T_E, f_trueRE.Clone("f_RE"), f_trueE.Clone("f_E"))

