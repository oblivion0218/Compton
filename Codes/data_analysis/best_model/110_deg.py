import ROOT
import numpy as np
import matplotlib.pyplot as plt
from lib import MoraPyRoot as mpr
from lib import LabLibrary as ll
import fit_models as fm


data_file_path = "/mnt/c/Users/User/Desktop/info/Compton/Measurments/Measurments_riflection/110_deg/"

H = ll.create_hist(data_file_path, "hist_sum.png")
Compton_peak = ll.search_photopeak(H, 0.4, 2)
sigma = 50

step = 30
n_steps = 3

min_fit = Compton_peak - (n_steps + 1) * step
max_fit = Compton_peak + n_steps * step

results_file_path = "/mnt/c/Users/User/Desktop/info/Compton/Measurments/Measurments_riflection/110_deg/plots/fit/best_model/"

fit_result_exp, f_background_exp, f_true_exp = fm.fit_exp_back(H, 100, Compton_peak, sigma, min_fit, max_fit, "Energy [keV]", "Counts", results_file_path)
fit_result_p1, f_background_p1, f_true_p1 = fm.fit_p1_back(H, 100, Compton_peak, sigma, min_fit, max_fit, "Energy [keV]", "Counts", results_file_path)
fit_result_p2, f_background_p2, f_true_p2 = fm.fit_p2_back(H, 100, Compton_peak, sigma, min_fit, max_fit, "Energy [keV]", "Counts", results_file_path)
fit_result_p4, f_background_p4, f_true_p4 = fm.fit_p4_back(H, 100, Compton_peak, sigma, min_fit, max_fit, "Energy [keV]", "Counts", results_file_path)



