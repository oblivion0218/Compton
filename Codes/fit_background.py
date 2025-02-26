import ROOT
import numpy as np
from lib import MoraPyRoot as mpr
from lib import LabLibrary as ll
   
data = ll.read_histogram_data("misura_22Na_20min_fondo.Spe")
hist = ROOT.TH1D(f"h", f"h", len(data), 0, len(data))

for i, bin_value in enumerate(data):
    hist.SetBinContent(i + 1, bin_value)

integral = hist.Integral()
if integral != 0:
    hist.Scale(1.0 / integral)
    
mpr.plot_hist_MPL(hist, "hist.png")

f_true = ROOT.TF1("f", "expo(0)", 0, 2000)
f_true.SetParameter(0, 8.4)
f_true.SetParameter(1, -0.012)

coo = [0.5, 0.9, 0.9, 0.6]
str0 = ["f1", "f2"]

mpr.stampa_graph_fit(hist, f_true, "fit2.png", "", "", "Counts", "", 55, 2000, 2, coo, str0)
mpr.plot_TF1_MPL(f_true, 0, 2000, "f.png")