import numpy as np
import matplotlib.pyplot as plt
import ROOT
from ROOT import TF1, TH1F, TCanvas, TLegend

# Lettura del file
def read_data(filename):
    data = []
    with open(filename, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if line.strip() and not line.startswith("Canale"):  # Ignora le prime righe non numeriche
                try:
                    _, _, eventi = map(int, line.split())
                    data.append(eventi)
                except ValueError:
                    continue
    return np.array(data)

# Creazione dell'istogramma ROOT
def create_histogram(data):
    hist = TH1F("hist", "Istogramma dei dati", len(data), 0, len(data))
    for i, value in enumerate(data):
        hist.SetBinContent(i+1, value)  # ROOT usa bin da 1
    return hist

# Funzione di fit: doppia gaussiana + fondo parabolico
def double_gauss_parabola(x, par):
    g1 = par[0] * np.exp(-0.5 * ((x[0] - par[1]) / par[2])**2)
    g2 = par[3] * np.exp(-0.5 * ((x[0] - par[4]) / par[5])**2)
    bg = par[6] + par[7]*x[0] + par[8]*x[0]**2
    return g1 + g2 + bg

# Creazione della funzione di fit
fit_func = TF1("fit_func", double_gauss_parabola, 1000, 1900, 9)

def setup_fit_function(fit_func , i):
    fit_func.SetParameters(150, 145, 100, 200, 1500, 100)
    fit_func.SetParLimits(0, 100, 200)
    fit_func.SetParLimits(1, 1400, 1500)
    fit_func.SetParLimits(2, 10, 120)
    fit_func.SetParLimits(3, 50,200)   #DEVO DEFINIRE MEGLIO QUESTI LIMITI!!!
    fit_func.FixParameter(4, 1490+i)
    fit_func.SetParLimits(5, 10,200)
    return fit_func

filename = "tot_15.txt"
data = read_data(filename)
hist = create_histogram(data)
hist.GetXaxis().SetRangeUser(1100,1800)

centro = []
err_centro = []

for i in range(20) :

    fit_func = setup_fit_function(fit_func , i )
    hist.Fit(fit_func, "R")



    #Creazione del canvas e fit
    canvas = TCanvas("c1", f"Fit Istogramma {1490 + i}" , 800, 600)
    hist.Draw()

    # Creazione delle singole componenti per il disegno
    gauss1 = TF1("gauss1", "[0]*exp(-0.5*((x-[1])/[2])**2)", 700, 2000)
    gauss1.SetParameters(fit_func.GetParameter(0), fit_func.GetParameter(1) , fit_func.GetParameter(2))
    gauss1.SetLineColor(ROOT.kGreen)

    gauss2 = TF1("gauss2", "[0]*exp(-0.5*((x-[1])/[2])**2)", 700, 2000)
    gauss2.SetParameters(fit_func.GetParameter(3), fit_func.GetParameter(4) , fit_func.GetParameter(5))
    gauss2.SetLineColor(ROOT.kBlue)

    bkg = TF1("bkg", "[0] + [1]*x + [2]*x*x", 700, 2000)
    bkg.SetParameters(fit_func.GetParameter(6), fit_func.GetParameter(7), fit_func.GetParameter(8))
    bkg.SetLineColor(ROOT.kOrange)
    bkg.SetLineStyle(2)  # Linea tratteggiata


    # Disegna le componenti individuali
    gauss1.Draw("SAME")
    gauss2.Draw("SAME")
    bkg.Draw("SAME")

    # Legenda
    legend = TLegend(0.6, 0.7, 0.9, 0.9)
    legend.AddEntry(hist, "Dati", "l")
    legend.AddEntry(fit_func, "tot", "l")
    legend.AddEntry(gauss1, "Compton 15°", "l")
    legend.AddEntry(gauss2, " 511 keV", "l")
    legend.AddEntry(bkg, "Fondo", "l")
    legend.Draw()

    canvas.Update()
    canvas.SaveAs(f"Fit Istogramma {1490 + i}.png")

  
    
    centro.append(fit_func.GetParameter(1))
    err_centro.append(fit_func.GetParError(1))
    

plt.hist(centro, bins=30, color='red', alpha=0.6, label='distribuzione centroidi', edgecolor='black')
plt.xlabel("centroidi")
plt.ylabel("Frequenza")
plt.legend()
plt.savefig("histo.png")
plt.show()