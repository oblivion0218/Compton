import os
import ROOT
from . import MoraPyRoot as mpr

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
        canvas.SaveAs(fileName + "_spectrum" + ".png")

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
        hist = ROOT.TH1D(f"h{j}", f"h{j}", len(data), 0, len(data))
        
        for i, bin_value in enumerate(data):
            hist.SetBinContent(i + 1, bin_value)
        
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
