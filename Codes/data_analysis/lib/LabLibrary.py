import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os
import ROOT
from . import MoraPyRoot as mpr

#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
# File Operations
#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
def get_file_names(directory_path):
    """
    Retrieves the names of all files in the specified directory.
    
    :param directory_path: Path to the directory containing the files.
    :return: List of file names in the directory.
    """
    file_names = []
    
    # Check if the provided path is a directory
    if os.path.isdir(directory_path):
        # Iterate through the items in the directory
        for file_name in os.listdir(directory_path):
            # Create the full file path
            full_file_path = os.path.join(directory_path, file_name)
            # Add the file name if it is a file (not a directory)
            if os.path.isfile(full_file_path):
                file_names.append(file_name)
    else:
        print(f"Il percorso {directory_path} non è una cartella valida.")
    
    return file_names


def remove_extension(file_list):
    """
    Removes file extensions from a list of file names.
    
    :param file_list: List of file names with extensions.
    :return: List of file names without extensions.
    """
    return [os.path.splitext(file)[0] for file in file_list]


#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
# Peak Search and Identification
#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
def search_peak(hist, noise_threshold, n_peaks, fileName=None):
    """
    Searches for peaks in a histogram using the TSpectrum algorithm.
    
    :param hist: ROOT histogram object.
    :param noise_threshold: Minimum threshold to identify peaks.
    :param n_peaks: Expected number of peaks to be found.
    :param fileName: (Optional) If provided, saves a plot of the histogram with peaks.
    :return: List of peak positions found in the histogram.
    """
    spectrum = ROOT.TSpectrum()
    n_found_peaks = spectrum.Search(hist, n_peaks, "", noise_threshold)

    if n_found_peaks == 0:
        raise Exception("Nessun picco trovato")
    
    # Get the positions of the peaks
    peak_positions = spectrum.GetPositionX()
    peak_positions_list = [peak_positions[i] for i in range(n_found_peaks)]

    print(f"Numero di picchi trovati: {n_found_peaks}")
    for i in range(n_found_peaks):
        print(f"Picco {i + 1}: posizione = {peak_positions[i]}")

    if fileName != None:
        canvas = ROOT.TCanvas("c1", "Istogramma con Picchi", 800, 600)
        hist.Draw() 
        canvas.SaveAs(fileName)

    return peak_positions_list


def search_photopeak(hist, noise_threshold, n_peaks, fileName=None):
    """
    Identifies the photopeak in a histogram, typically for 511 keV photons from Na-22.
    
    :param hist: ROOT histogram object.
    :param noise_threshold: Minimum threshold to identify peaks.
    :param n_peaks: Expected number of peaks to be found.
    :param fileName: (Optional) If provided, saves a plot of the histogram with peaks.
    :return: Position of the photopeak in the histogram.
    """
    peak_positions = search_peak(hist, noise_threshold, n_peaks, fileName)
    max_position = peak_positions[0]
    
    n_found_peaks = len(peak_positions)
    # Find the maximum peak position
    for i in range(n_found_peaks):
        print(f"Picco {i + 1}: posizione = {peak_positions[i]}")
        if peak_positions[i] > max_position:
            max_position = peak_positions[i]
    
    return max_position


def search_first_peak(hist, noise_threshold, n_peaks, fileName=None):  
    """
    Identifies the first peak in a histogram.
    
    :param hist: ROOT histogram object.
    :param noise_threshold: Minimum threshold to identify peaks.
    :param n_peaks: Expected number of peaks to be found.
    :param fileName: (Optional) If provided, saves a plot of the histogram with peaks.
    :return: Position of the first peak in the histogram.
    """
    peak_positions = search_peak(hist, noise_threshold, n_peaks, fileName)
    min_position = peak_positions[0]
    
    n_found_peaks = len(peak_positions)
    # Find the maximum peak position
    for i in range(n_found_peaks):
        print(f"Picco {i + 1}: posizione = {peak_positions[i]}")
        if peak_positions[i] < min_position:
            min_position = peak_positions[i]
    
    return min_position


#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
# Calibration
#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
def calibration(x, a = 0.36, b = -11.5):
    """
    Calibration function for converting ADC channels to energy in keV.
    
    :param x: ADC channel value.
    :param a: Calibration coefficient (default: 0.36).
    :param b: Calibration offset (default: -11.5).
    :return: Energy in keV.
    """
    return a * x + b


#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
# Fit functions for Compton study
#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
def fit_photopeak_linear_background(hist, fileNamePNG, noise_threshold, n_peaks):
    """
    Fits the photopeak in a histogram with a Gaussian function over a linear background.
    
    :param hist: ROOT histogram object.
    :param fileNamePNG: Path to save the plot of the fit.
    :param noise_threshold: Minimum threshold to identify peaks.
    :param n_peaks: Expected number of peaks to be found.
    :return: Tuple containing the fitted histogram and the fit function.
    """
    coo2 = [0.6, 0.13, 0.9, 0.53]
    str2 = ["f1", "f2", "Amp", "<x>", "#sigma"]

    photopeak_x = search_photopeak(hist, noise_threshold, n_peaks)
    
    n_bins = hist.GetNbinsX()

    # Determina l'estremo per il fit
    if n_bins >= 1800: 
        extreme = hist.GetXaxis().GetBinUpEdge(n_bins) / 10
    else: 
        extreme = hist.GetXaxis().GetBinUpEdge(n_bins) / 50

    # Gaussian fit
    f_picco = ROOT.TF1("picco", "gaus(0)", 0, 2000)
    # L'inizializzazione dei parametri varia al variare del numero di bin
    f_picco.SetParameters(0, 40000 / n_bins)
    f_picco.SetParameters(1, photopeak_x)
    f_picco.SetParameters(2, 0.001 * n_bins)
    mpr.fit(hist, f_picco, 3, "S", 1000, photopeak_x - extreme, photopeak_x + extreme)

    # Linear fit
    f_fondo = ROOT.TF1("fondo", "pol1(0)", 0, 2000)
    f_fondo.SetParameters(0, 1)
    f_fondo.SetParameters(1, 1)
    mpr.fit(hist, f_fondo, 3, "S", 1000, photopeak_x - extreme * 2, photopeak_x + extreme * 2)

    # Background + peak
    f_true = ROOT.TF1("modello", "pol1(0) + gaus(2)", 0, 2000)
    f_true.SetParameter(0, f_fondo.GetParameter(0))
    f_true.SetParameter(1, f_fondo.GetParameter(1))
    f_true.SetParameter(2, f_picco.GetParameter(0))
    f_true.SetParameter(3, f_picco.GetParameter(1))
    f_true.SetParameter(4, f_picco.GetParameter(2))

    extreme_graph = ("x", [photopeak_x - extreme, photopeak_x + extreme])

    mpr.stampa_graph_fit_range(hist, f_true, extreme_graph, fileNamePNG, "", "", "Counts", "", photopeak_x - extreme, photopeak_x + extreme, 5, coo2, str2)


# For histograms normalized to have Area = 1
def stampa_graph_fit_ComptonStudy(hist, f_true, scale_factor, min, max, file_path, fileNamePNG, graph_name, x_axis_name, y_axis_name, graphic_option, pave_coordinates = None, f_background = None):
    """
    Fit a peak with a given function and print the number of events, the mean energy, and the energy resolution.
    
    :param hist: ROOT histogram object.
    :param f_true: Fit function.
    :param scale_factor: Scale factor for the number of events.
    :param min: Minimum value for the fit.
    :param max: Maximum value for the fit.
    :param file_path: Path to save the plot.
    :param fileNamePNG: Name of the output file.
    :param graph_name: Name of the graph.
    :param x_axis_name: Name of the x-axis.
    :param y_axis_name: Name of the y-axis.
    :param graphic_option: Graphic option.
    :param pave_coordinates: Coordinates of the text box.
    :param f_background: Background function.
    """
    canvas = ROOT.TCanvas()
    
    hist.Draw(graphic_option)
    hist.SetTitle(graph_name)
    hist.GetXaxis().SetTitle(x_axis_name)
    hist.GetYaxis().SetTitle(y_axis_name)

    fit_result = hist.Fit(f_true, "S", "", min, max) 

    if pave_coordinates != None: 
        text_box = ROOT.TPaveText(pave_coordinates[0], pave_coordinates[1], pave_coordinates[2], pave_coordinates[3], "NDC")
        text_box.SetFillColor(0)
        text_box.SetTextAlign(12)

        E_mean = f_true.GetParameter(1)
        E_mean_error = f_true.GetParError(1)
        E_mean_SN = mpr.exponential(E_mean_error)

        if abs(E_mean_SN.exp) < 3:
            text_box.AddText(f"<E> = {E_mean:.3f} #pm {E_mean_error:.3f}")
        else:   
            text_box.AddText(f"E_mean = ({E_mean * 10 ** -E_mean_SN.exp:.3f} #pm {E_mean_SN.n:.3f}) * 10^{{{E_mean_SN.exp}}}")

        FWHM = 2.355 * f_true.GetParameter(2)

        ER = FWHM / E_mean
        ER_error = np.sqrt((2.355 * f_true.GetParError(2) / E_mean)**2 + (FWHM * f_true.GetParError(1) / E_mean**2)**2)
        ER_SN = mpr.exponential(ER_error)

        if abs(ER_SN.exp) < 3:
            text_box.AddText(f"ER = {ER:.3f} #pm {ER_error:.3f}")
        else:    
            text_box.AddText(f"ER = ({ER * 10 ** -ER_SN.exp:.3f} #pm {ER_SN.n:.3f}) * 10^{{{ER_SN.exp}}}")

        if f_background:
            N_hit_pc = f_true.Integral(min, max) - f_background.Integral(min, max)
        else:
            N_hit_pc = f_true.Integral(min, max)

        N_hit = N_hit_pc * scale_factor
        N_hit_error = np.sqrt(N_hit)
        N_hit_SN = mpr.exponential(N_hit_error)

        if abs(N_hit_SN.exp) < 3:
            text_box.AddText(f"N = {N_hit:.3f} #pm {N_hit_error:.3f}")
        else:
            text_box.AddText(f"N = ({N_hit * 10 ** -N_hit_SN.exp:.3f} #pm {N_hit_SN.n:.3f}) * 10^{{{N_hit_SN.exp}}}")

        N_hit_pc_error = N_hit_error / scale_factor
        N_hit_pc_SN = mpr.exponential(N_hit_pc_error)

        if abs(N_hit_pc_SN.exp) < 3:
            text_box.AddText(f"N% = {N_hit_pc:.3f} #pm {N_hit_pc_error:.3f}")
        else:
            text_box.AddText(f"N% = ({N_hit_pc * 10 ** -N_hit_pc_SN.exp:.3f} #pm {N_hit_pc_SN.n:.3f}) * 10^{{{N_hit_pc_SN.exp}}}")

        text_box.Draw()

        canvas.Print(file_path + "plots/fit/" + fileNamePNG, "png")
        del text_box
    
    else:            
        canvas.Print(file_path + "plots/fit/" + fileNamePNG, "png")

def compton_scattering(theta):
    return 511 / (2 - np.cos(theta))

def efficiency (Angle):

    A = 1.55710
    B = -0.09831
    C = 3.53226
    D = 0.10209

    Energy = compton_scattering(Angle)

    return  A* Energy^(-B) * np.exp(-Energy*C) + D

def plot_results(hist, fit_result, f_background, f_true, rebin_param, min_fit, max_fit, file_path, fileNamePNG, x_axis_name, y_axis_name, time):
    """
    Plot the results of the fit."

    :param hist: ROOT histogram object."
    :param fit_result: Fit result."
    :param f_background: Background function."
    :param f_true: True function."
    :param rebin_param: Rebinning parameter."
    :param min_fit: Minimum value of the fit."
    :param max_fit: Maximum value of the fit."
    :param file_path: Path where to save the plot."
    :param fileNamePNG: Name of the graph."
    :param x_axis_name: Name of the x-axis."
    :param y_axis_name: Name of the y-axis."
    """
    chi2 = fit_result.Chi2()
    ndf = fit_result.Ndf()
    chi2_ndf = chi2 / ndf

    E_mean = (f_true.GetParameter(3), f_true.GetParError(3))
    sigma = (f_true.GetParameter(4), f_true.GetParError(4))
    FWHM = (2.355 * sigma[0], 2.355 * sigma[1])
    ER = (FWHM[0] / E_mean[0], np.sqrt((FWHM[1] / E_mean[0]) ** 2 + (FWHM[0] * E_mean[1] / E_mean[0] ** 2) ** 2))

    integral = (f_true.Integral(min_fit, max_fit) - f_background.Integral(min_fit, max_fit)) / rebin_param
    N_hit = (integral, np.sqrt(integral))
    rate = (N_hit[0] / time, N_hit[1] / time)

    text = rf"$\chi^{{2}}/\mathrm{{dof}} = {chi2:.3f}/{ndf} = {chi2_ndf:.3f}$"
    text += "\n"
    text += f"<E> = {E_mean[0]:.2f} ± {E_mean[1]:.2f}\n"
    text += f"ER = {ER[0]:.3f} ± {ER[1]:.3f}\n"
    text += f"N = {N_hit[0]:.2f} ± {N_hit[1]:.2f}\n"
    text += f"Rate = {rate[0]:.4f} ± {rate[1]:.4f} Hz\n"

    #-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    # PLOT RESULTS - total fit and background
    #-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-       
    n_bins = hist.GetNbinsX()
    bin_edges = np.array([hist.GetBinLowEdge(i+1) for i in range(n_bins+1)])
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    bin_values = np.array([hist.GetBinContent(i+1) for i in range(n_bins)])

    x_values = np.linspace(min_fit, max_fit, 1000)
    f_true_func = np.vectorize(lambda x: f_true.Eval(x))
    f_back_func = np.vectorize(lambda x: f_background.Eval(x))

    y_true = f_true_func(x_values)
    y_back = f_back_func(x_values)
    
    # Create a figure with two subplots: one for the main plot and one for the residuals
    fig = plt.figure(figsize=(8, 10))
    gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1.5])  
    
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
        if min_fit <= center <= max_fit:
            model_value = f_true.Eval(center)
            residual = (bin_values[i] - model_value) / model_value
            residuals.append(residual)
            residual_centers.append(center)  # Only include centers in the integration domain

    ax_residual.errorbar(residual_centers, residuals, fmt='x', color='black', label="Residuals")
    ax_residual.axhline(0, color="red", linewidth=2, label="Zero Line")

    ax_residual.set_xlim(min_fit - 200, max_fit + 200)
    ax_residual.set_ylabel(r"$\frac{data - model}{model}$", size=12)
    ax_residual.set_xlabel(x_axis_name)
    ax_residual.grid(True, which='both', linestyle='--', linewidth=0.5)
    ax_residual.legend(loc="upper right")

    plt.tight_layout()
    plt.savefig(file_path + fileNamePNG)
    plt.close()

    return N_hit, rate


#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
# Histogram Operations
#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
def read_histogram_data(filename):
    """
    Reads histogram data from a .Spe file.
    
    :param filename: Path to the .Spe file.
    :return: List of bin values from the histogram.
    """
    with open(filename, "r") as file:
        histogram_data = []
        found_data_section = False
        
        # Search for the $DATA section
        for line in file:
            if "$DATA:" in line:
                found_data_section = True
                break
        
        if not found_data_section:
            raise Exception("Section $DATA not found in file")
        
        # Skip the line with the range "0 2047"
        next(file)
        
        # Read the histogram bin values
        for line in file:
            try:
                bin_value = int(line.strip())
                histogram_data.append(bin_value)
            except ValueError:
                break

    return histogram_data


def hist_vector(directory_path):
    """
    Reads multiple histogram files from a directory and returns a list of ROOT histograms.
    
    :param directory_path: Path to the directory containing histogram files.
    :return: List of ROOT histogram objects.
    """
    hist_file_names = get_file_names(directory_path)
    hist_list = []

    for j, hist_file_name in enumerate(hist_file_names):
        
        data = read_histogram_data(directory_path + hist_file_name)
        n_bins = int(len(data))
        hist = ROOT.TH1D(f"h{j}", f"h{j}", n_bins, 0, n_bins)
        
        for i, bin_value in enumerate(data):
            hist.SetBinContent(int(i) + 1, bin_value)
        
        hist_list.append(hist)
    
    return hist_list


def spectum_sum(histograms):
    """
    Computes the sum of multiple histograms.
    
    :param histograms: List of ROOT histogram objects.
    :return: ROOT histogram representing the sum of all input histograms.
    """
    hist_sum = ROOT.TH1F("hist_sum", "Sum of spectra", histograms[0].GetNbinsX(), histograms[0].GetXaxis().GetXmin(), histograms[0].GetXaxis().GetXmax())
    
    for hist in histograms:
        hist_sum.Add(hist)
    
    return hist_sum


def create_hist(file_path, fileNamePNG):
    hs = hist_vector(file_path + "data/")

    for h in hs:
        mpr.plot_hist_MPL(h, file_path + "plots/hist/" + h.GetName() + ".png")
    
    H = spectum_sum(hs)
    mpr.plot_hist_MPL(H, file_path + "plots/hist/" + fileNamePNG)

    return H


def normalize_histogram(hist, normalized_area=1):
    """
    Normalizes a histogram to have an area of 1.
    
    :param hist: ROOT histogram object.
    :return: Normalized ROOT histogram object.
    """
    integral = hist.Integral()
    scale_factor = normalized_area / integral
    for i in range(1, hist.GetNbinsX() + 1):
        hist.SetBinContent(i, hist.GetBinContent(i) * scale_factor)
    
    return hist


#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
# Fit stability study
#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
def stability_study_extreme(fit_peaks, H, peakCompton, sigmaCompton, step, max_step, x_axis_name, y_axis_name, file_path):
    """
    Study the stability of the fit as a function of the fitting range extremes.
    In addition to χ²/ndf, also evaluate the fitted Compton peak position, σ, and
    the energy resolution (ER = σ/peak × 100%).
    
    :param fit_peaks: Function to fit the peaks, with arguments: 
            :param H: ROOT histogram object.
            :param peakCompton: Nominal peak position.
            :param sigmaCompton: Nominal sigma of the peak.
            :param min_fit: Minimum value of the fit.
            :param max_fit: Maximum value of the fit.
            :param x_axis_name: Label for x-axis.
            :param y_axis_name: Label for y-axis.
            :param file_path: Path (folder) to save the plots.
    :param H: ROOT histogram object.
    :param peakCompton: Nominal peak position.
    :param sigmaCompton: Nominal sigma of the peak.
    :param step: Step size for changing the fit extremes.
    :param max_step: Maximum step (number of iterations).
    :param x_axis_name: Label for x-axis.
    :param y_axis_name: Label for y-axis.
    :param file_path: Path (folder) to save the plots.
    """
    n_steps = [i for i in range(1, max_step)]
    chi2_list = []
    peak_fit_list = []
    sigma_fit_list = []
    er_list = []
    
    # Loop over fit extremes
    for i in n_steps:
        min_fit = peakCompton - i * step
        max_fit = peakCompton + i * step
        # Call your fit_peaks routine (which returns a canvas, background fit, and the full fit function)
        c, _, f_true = fit_peaks(H, peakCompton, sigmaCompton,
                                 min_fit, max_fit,
                                 x_axis_name, y_axis_name,
                                 file_path + "fit_stability/extreme/" + str(i))
        chi2_val = c.Chi2() / c.Ndf()
        fitted_peak = f_true.GetParameter(3)  # Compton peak position
        fitted_sigma = f_true.GetParameter(4)  # Compton sigma
        er = (fitted_sigma * 2.355 / fitted_peak) * 100  # Energy resolution in %
        
        chi2_list.append(chi2_val)
        peak_fit_list.append(fitted_peak)
        sigma_fit_list.append(fitted_sigma)
        er_list.append(er)
    
    # Plot all the stability studies in a single figure using subplots.
    fig, axs = plt.subplots(2, 2, figsize=(18, 10))
    
    # Subplot 1: χ²/ndf stability
    axs[0, 0].plot(n_steps, chi2_list, marker='o', linestyle='--', color='blue')
    axs[0, 0].set_title(r'$\chi^2/ndf$ Stability')
    axs[0, 0].set_xlabel("N_step")
    axs[0, 0].set_ylabel(r"$\chi^2/ndf$")
    axs[0, 0].grid(True)
    axs[0, 0].set_xticks(n_steps)
    
    # Subplot 2: Fitted Peak Position
    axs[0, 1].plot(n_steps, peak_fit_list, marker='o', linestyle='--', color='red')
    axs[0, 1].set_title("Fitted Compton Peak Position")
    axs[0, 1].set_xlabel("N_step")
    axs[0, 1].set_ylabel("Peak Position")
    axs[0, 1].grid(True)
    axs[0, 1].set_xticks(n_steps)
    
    # Subplot 3: Fitted Sigma
    axs[1, 0].plot(n_steps, sigma_fit_list, marker='o', linestyle='--', color='green')
    axs[1, 0].set_title("Fitted σ Stability")
    axs[1, 0].set_xlabel("N_step")
    axs[1, 0].set_ylabel("σ")
    axs[1, 0].grid(True)
    axs[1, 0].set_xticks(n_steps)
    
    # Subplot 4: Energy Resolution (ER = σ/peak × 100)
    axs[1, 1].plot(n_steps, er_list, marker='o', linestyle='--', color='magenta')
    axs[1, 1].set_title("Energy Resolution Stability")
    axs[1, 1].set_xlabel("N_step")
    axs[1, 1].set_ylabel("ER (%)")
    axs[1, 1].grid(True)
    axs[1, 1].set_xticks(n_steps)
    
    fig.tight_layout()
    plot_file = file_path + "stability_summary_extreme.png"
    plt.savefig(plot_file)
    plt.close()
    print(f"Stability summary plot saved to {plot_file}")


def stability_study_rebin(fit_peaks, H, peakCompton, sigmaCompton, rebin_max, min_fit, max_fit, x_axis_name, y_axis_name, file_path):
    """
    Study the stability of the fit as a function of the histogram rebinning factor.
    In addition to χ²/ndf, evaluate the fitted Compton peak position, σ, and the energy 
    resolution (ER = σ/peak × 100%).
    
    :param fit_peaks: Function to fit the peaks, with arguments: 
            :param H: ROOT histogram object.
            :param peakCompton: Nominal peak position.
            :param sigmaCompton: Nominal sigma of the peak.
            :param min_fit: Minimum value of the fit.
            :param max_fit: Maximum value of the fit.
            :param x_axis_name: Label for x-axis.
            :param y_axis_name: Label for y-axis.
            :param file_path: Path (folder) to save the plots.
    :param H: ROOT histogram object.
    :param peakCompton: Nominal peak position.
    :param sigmaCompton: Nominal sigma of the peak.
    :param rebin_max: Maximum rebin factor to loop over (from 1 to rebin_max).
    :param min_fit: Minimum value of the fit.
    :param max_fit: Maximum value of the fit.
    :param x_axis_name: Label for x-axis.
    :param y_axis_name: Label for y-axis.
    :param file_path: Path (folder) to save the plots.
    """

    # Create lists to store stability parameters
    rebin_factors = [i for i in range(1, rebin_max + 1)]
    chi2_list = []
    peak_fit_list = []
    sigma_fit_list = []
    er_list = []
    
    # Loop over rebin factors
    for i in rebin_factors:
        # Clone and rebin the histogram to avoid modifying the original
        H_rebin = H.Clone()
        H_rebin.Rebin(i)
        
        # Perform the fit on the rebinned histogram
        c, _, f_true = fit_peaks(H_rebin, peakCompton, sigmaCompton,
                                 min_fit, max_fit,
                                 x_axis_name, y_axis_name,
                                 file_path + "fit_stability/rebin/" + str(i))
        # Extract χ²/ndf
        chi2_val = c.Chi2() / c.Ndf()
        # Extract fitted parameters (using indices as in your fit_peaks)
        fitted_peak = f_true.GetParameter(3)  # Compton peak position
        fitted_sigma = f_true.GetParameter(4)   # Compton sigma
        er = (fitted_sigma * 2.355 / fitted_peak) * 100   # Energy resolution in %
        
        chi2_list.append(chi2_val)
        peak_fit_list.append(fitted_peak)
        sigma_fit_list.append(fitted_sigma)
        er_list.append(er)
    
    # Plot all the stability studies in a single figure using subplots
    fig, axs = plt.subplots(2, 2, figsize=(18, 10))
    
    # Subplot 1: χ²/ndf stability vs. rebin factor
    axs[0, 0].plot(rebin_factors, chi2_list, marker='o', linestyle='--', color='blue')
    axs[0, 0].set_title(r'$\chi^2/ndf$ Stability vs Rebin Factor')
    axs[0, 0].set_xlabel("Rebin Factor")
    axs[0, 0].set_ylabel(r"$\chi^2/ndf$")
    axs[0, 0].grid(True)
    axs[0, 0].set_xticks(rebin_factors)
    
    # Subplot 2: Fitted Peak Position vs. rebin factor
    axs[0, 1].plot(rebin_factors, peak_fit_list, marker='o', linestyle='--', color='red')
    axs[0, 1].set_title("Fitted Compton Peak Position vs Rebin Factor")
    axs[0, 1].set_xlabel("Rebin Factor")
    axs[0, 1].set_ylabel("Peak Position")
    axs[0, 1].grid(True)
    axs[0, 1].set_xticks(rebin_factors)
    
    # Subplot 3: Fitted σ vs. rebin factor
    axs[1, 0].plot(rebin_factors, sigma_fit_list, marker='o', linestyle='--', color='green')
    axs[1, 0].set_title("Fitted σ Stability vs Rebin Factor")
    axs[1, 0].set_xlabel("Rebin Factor")
    axs[1, 0].set_ylabel("σ")
    axs[1, 0].grid(True)
    axs[1, 0].set_xticks(rebin_factors)
    
    # Subplot 4: Energy Resolution vs. rebin factor
    axs[1, 1].plot(rebin_factors, er_list, marker='o', linestyle='--', color='magenta')
    axs[1, 1].set_title("Energy Resolution vs Rebin Factor")
    axs[1, 1].set_xlabel("Rebin Factor")
    axs[1, 1].set_ylabel("ER (%)")
    axs[1, 1].grid(True)
    axs[1, 1].set_xticks(rebin_factors)
    
    fig.tight_layout()
    plot_file = file_path + "stability_summary_rebin.png"
    plt.savefig(plot_file)
    plt.close()
    print(f"Rebin stability summary plot saved to {plot_file}")
