import ROOT
import glob
import os

# Funzione per fittare un picco con una gaussiana e fondo lineare
def fit_gaussian_with_linear_background(histogram, peak_position, output_file, hist_name):
    # Definisci la funzione: gaussiana + fondo lineare
    func = ROOT.TF1("fit_func", "[0] + [1]*x + [2]*exp(-0.5*((x - [3])/[4])**2)", peak_position-200, peak_position+200)
    
    # Parametri iniziali per il fit
    func.SetParameters(1, -0.1, 220, peak_position, 20)  # [0]=offset lineare, [1]=pendenza, [2]=ampiezza, [3]=centro gaussiana, [4]=sigma
    # Vincolo sulla posizione del picco (centro della gaussiana) tra 500 e 800
    func.SetParLimits(3, 500, 800)
    # Esegui il fit
    histogram.Fit(func, "R")
    # Disegna la funzione del fit sovrapposta all'istogramma
    func.Draw("same")
    
    # Stampa i risultati del fit nel file di output
    output_file.write(f"{hist_name}")
    for i in range(5):
        param = func.GetParameter(i)
        error = func.GetParError(i)
        output_file.write(f"\t{param:.4f}\t{error:.4f}")
    output_file.write("\n")

#----------------MAIN----------------------

# Directory dove si trovano i file
directory_path = ""
# Pattern per trovare i file, ad esempio per trovare i file che iniziano con un numero
file_pattern = os.path.join(directory_path, "[0-9]*.Spe")

# Apri il file di output per scrivere i risultati del fit
with open("fit_results.txt", "w") as output_file:
    # Scrivi l'intestazione del file
    output_file.write("hist_name\tOffset\tOffset_err\tSlope\tSlope_err\tAmplitude\tAmplitude_err\tCenter\tCenter_err\tSigma\tSigma_err\n")
    
    # Itera su tutti i file che corrispondono al pattern
    for file_path in glob.glob(file_pattern):
        # Lista per contenere i valori
        values = []

        # Leggi il file saltando le prime 12 righe e le ultime 15
        with open(file_path, 'r') as file:
            lines = file.readlines()[12:-15]  

        # Itera sulle righe rimanenti
        for line in lines:
            # Rimuovi eventuali spazi bianchi e salta righe vuote
            line = line.strip()
            if line:
                # Converte il valore in float e aggiungilo all'array
                values.append(float(line))

        # Crea un istogramma per il file corrente
        if values:  # Assicura che ci siano dati validi
            hist_name = os.path.basename(file_path).replace(".Spe", "")
            histogram = ROOT.TH1F(hist_name + "UGO", f"Istogramma di {hist_name}", 2048, 0, 2048)
            for bin_index, value in enumerate(values):
                histogram.SetBinContent(bin_index + 1, value)  # `SetBinContent` accetta l'indice del bin a partire da 1
            
            # Crea un canvas per visualizzare l'istogramma
            canvas = ROOT.TCanvas(f"canvas_{hist_name}", f"Istogramma {hist_name} UGO", 800, 600)
            # Imposta il range dell'asse X tra 500 e 800
            histogram.GetXaxis().SetRangeUser(500, 800)
            histogram.Draw("HIST")  # Usa l'opzione "HIST" per il riempimento completo
            
            # Determina la posizione del picco in base al nome dell'istogramma
            if hist_name == "0,5":
                peak_position = 570
            else:
                peak_position = 650

            # Esegui il fit con la funzione gaussiana più fondo lineare
            fit_gaussian_with_linear_background(histogram, peak_position, output_file, hist_name)

            # Salva il canvas come immagine PNG per riferimento
            canvas.SaveAs(f"SHAPING{hist_name}UGO.png")
