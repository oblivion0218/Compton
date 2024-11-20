import os
import ROOT

def get_file_names(directory_path):
    # Lista per contenere i nomi dei file
    file_names = []
    
    # Controlla se il percorso fornito è una directory
    if os.path.isdir(directory_path):
        # Itera attraverso gli elementi della cartella
        for file_name in os.listdir(directory_path):
            # Crea il percorso completo del file
            full_file_path = os.path.join(directory_path, file_name)
            # Aggiunge il nome del file se è un file (non una cartella)
            if os.path.isfile(full_file_path):
                file_names.append(file_name)
    else:
        print(f"Il percorso {directory_path} non è una cartella valida.")
    
    return file_names

def remove_extension(file_list):
    # Usa os.path.splitext per separare il nome del file dall'estensione
    return [os.path.splitext(file)[0] for file in file_list]

# Funzione per leggere i dati dell'istogramma da file .Spe
def read_histogram_data(filename):
    with open(filename, "r") as file:
        histogram_data = []
        found_data_section = False
        
        # Cerca la sezione $DATA
        for line in file:
            if "$DATA:" in line:
                found_data_section = True
                break
        
        if not found_data_section:
            raise Exception("Section $DATA not found in file")
        
        # Salta la riga con l'intervallo "0 2047"
        next(file)
        
        # Legge i valori dei bin dell'istogramma
        for line in file:
            try:
                bin_value = int(line.strip())
                histogram_data.append(bin_value)
            except ValueError:
                break

    return histogram_data

# ritorna la posizione del photopeack for 511 keV photon in Na22 
def search_photopeak(hist, noise_threshold, n_peaks, fileName=None):
    spectrum = ROOT.TSpectrum()
    n_found_peaks = spectrum.Search(hist, n_peaks, "", noise_threshold)
    
    print(f"Numero di picchi trovati: {n_found_peaks}")
    
    if n_found_peaks == 0:
        raise Exception("Nessun picco trovato")
    
    # Ottieni le posizioni dei picchi trovati
    peak_positions = spectrum.GetPositionX()
    
    max_position = peak_positions[0]
    
    # Trova la posizione del picco massimo
    for i in range(n_found_peaks):
        print(f"Picco {i + 1}: posizione = {peak_positions[i]}")
        if peak_positions[i] > max_position:
            max_position = peak_positions[i]
    
    # Disegna l'istogramma e aggiungi marker sui picchi trovati
    if fileName != None:
        canvas = ROOT.TCanvas("c1", "Istogramma con Picchi", 800, 600)
        hist.Draw() 
        canvas.SaveAs(fileName + "_spectrum" + ".png")
    
    return max_position