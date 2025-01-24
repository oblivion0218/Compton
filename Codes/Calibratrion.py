import MoraPyRoot as mpr 
import LabLibrary as ll
import ROOT
import numpy as np

file_path = "/mnt/c/Users/User/Desktop/info/Compton/Calibration/"

def peak_error_N(peak):
    bin_peak = hist.FindBin(peak)
    bin_content = hist.GetBinContent(bin_peak)
    bin_width = hist.GetBinWidth(bin_peak)
    return bin_width / np.sqrt(bin_content)

def peak_error_FWHM(peak):
    bin_peak = hist.FindBin(peak)
    fwhm = hist.GetBinWidth(bin_peak) 
    return fwhm / 2.3548

# Ti 44
fileName = "44-Ti"
data = ll.read_histogram_data(file_path + "data/" + fileName + ".Spe")
hist = ROOT.TH1D("h", "h", len(data), 0, len(data))

for i in range(len(data)):
    hist.SetBinContent(i + 1, data[i])  # i+1 perché i bin in ROOT partono da 1

mpr.plot_hist_MPL(hist, file_path + "plots/hist/" +fileName + ".png")
peaks_positions = sorted(ll.search_peak(hist, noise_threshold=0.03, n_peaks=10, fileName=file_path + "plots/hist/" + fileName + "_peaks.png"))

for peak in peaks_positions:
    print(peak, peak_error_N(peak), peak_error_FWHM(peak))

print("\nFine\n")