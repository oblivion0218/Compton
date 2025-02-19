import ROOT
import numpy as np
from lib import MoraPyRoot as mpr
from lib import LabLibrary as ll

file_path = "/mnt/c/Users/User/Desktop/info/Compton/Measurments_trasmission/50_deg/"

hs = ll.hist_vector(file_path + "data/")

for h in hs:
    mpr.plot_hist_MPL(h, file_path + "plots/hist/" + h.GetName() + ".png")
H = ll.spectum_sum(hs)
mpr.plot_hist_MPL(H, file_path + "plots/hist/" + "hist_sum.png")

ll.search_photopeak(H, 0.1, 4, file_path + "plots/" + "511_peak.png")
ll.search_photopeak(H, 0.5, 4, file_path + "plots/" + "Compton_peak.png")
