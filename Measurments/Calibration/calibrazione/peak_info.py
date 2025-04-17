import ROOT
import LabLibrary as lab

def save_fit_results(filename, fit_results):
    """Salva i risultati del fit in un file."""
    with open(filename, 'w') as f:
        f.write("# Picco, Centro, Altezza \n" )
        for result in fit_results:
            f.write(", ".join(map(str, result)) + "\n")
            

def search_photopeak(hist, noise_threshold, n_peaks, fileName=None):
    spectrum = ROOT.TSpectrum()
    n_found_peaks = spectrum.Search(hist, n_peaks, "", noise_threshold)
    
    print(f"Numero di picchi trovati: {n_found_peaks}")
    
    # Ottieni le posizioni dei picchi trovati
    peak_positions = spectrum.GetPositionX()

    for i in range(n_found_peaks):
        print(f"Picco {i + 1}: posizione = {peak_positions[i]}")

    peak_positions_list = [peak_positions[i] for i in range(n_found_peaks)]
  
    # Disegna l'istogramma e aggiungi marker sui picchi trovati
    if fileName != None:
        canvas = ROOT.TCanvas("c1", "Istogramma con Picchi", 800, 600)
        hist.Draw() 
        canvas.SaveAs(fileName + "_spectrum" + ".png")
    
    return peak_positions_list


name = "44-Ti"
spe_file = name + ".Spe"
noise = 0.01
n_peaks = 10

print(f"Lettura dati da {spe_file}...")
histogram_data = lab.read_histogram_data(spe_file)

# Creazione dell'oggetto TH1F
n_bins = len(histogram_data)
hist = ROOT.TH1F("hist", "Spettro", n_bins, 0, n_bins)
for i, value in enumerate(histogram_data):
    hist.SetBinContent(i+1, value)

position = search_photopeak(hist, noise , n_peaks , name)

position.sort()

fit_results = []  # Lista per memorizzare i risultati del fit

for i in range(len(position)):

    amp_guess = hist.GetBinContent(hist.FindBin(position[i]))

    # Salvare i risultati in una lista
    fit_results.append([
        i + 1, position[i] , amp_guess
    ])


# Salvare i risultati del fit in un file di testo
save_fit_results(f"{name}_fit_results.txt", fit_results)


