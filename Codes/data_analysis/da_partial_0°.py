import ROOT
import numpy as np
from lib import MoraPyRoot as mpr
from lib import LabLibrary as ll

# === Percorso ai dati ===
file_path = "/mnt/c/Users/ASUS/Desktop/WSL_shared/Compton/Measurments/Measurments_riflection/0_deg/"
H = ll.create_hist(file_path, "hist_sum.png")  # Istogramma ROOT


# === Guess iniziali: (mu, sigma, ampiezza) ===
mu_guess = 1487     
sigma_guess = 55
amp_guess = 3500


# === Background parabola ===
a_guess, b_guess, c_guess = 0, 0, 600


# === Range di fit ===
fit_min = mu_guess - 3 * sigma_guess
fit_max = mu_guess + 3 * sigma_guess


# === Definizione della funzione di fit ===
fit_func = ROOT.TF1("fit_func", "[0]*exp(-0.5*((x-[1])/[2])^2) + [3]*x*x + [4]*x + [5]", fit_min, fit_max)
fit_func.SetParameters(amp_guess, mu_guess, sigma_guess, a_guess, b_guess, c_guess)


# === Fit ===
H.Fit(fit_func, "R")  # "R" = usa solo l'intervallo specificato


# === Parametri e loro errori ===
A = fit_func.GetParameter(0)
mu = fit_func.GetParameter(1)
sigma = fit_func.GetParameter(2)
sigma_A = fit_func.GetParError(0)
sigma_sigma = fit_func.GetParError(2)
a = fit_func.GetParameter(3)
b = fit_func.GetParameter(4)
c = fit_func.GetParameter(5)
sigma_a = fit_func.GetParError(3)
sigma_b = fit_func.GetParError(4)
sigma_c = fit_func.GetParError(5)


# === Calcolo area della gaussiana (sottesa) ===
area_totale = fit_func.Integral(fit_min, fit_max)
parabola_only = ROOT.TF1("parabola", "[0]*x*x + [1]*x + [2]", fit_min, fit_max)
parabola_only.SetParameters(a,b,c)
area_parabola = parabola_only.Integral(fit_min, fit_max)
area_net = area_totale - area_parabola


# === Calcolo dell'errore sull'area netta ===
factor = np.sqrt(2 * np.pi)
sigma_area_gauss = np.sqrt((sigma * factor * sigma_A)**2 + (A * factor * sigma_sigma)**2) 

# Intervallo di integrazione
x1 = fit_min
x2 = fit_max

# Derivate parziali
da = (x2**3 - x1**3) / 3
db = (x2**2 - x1**2) / 2
dc = (x2 - x1)

# Errore sul fondo
sigma_area_parabola = np.sqrt((da * sigma_a)**2 + (db * sigma_b)**2 + (dc * sigma_c)**2)

# Errore totale su area_net
sigma_area_net = np.sqrt(sigma_area_gauss**2 + sigma_area_parabola**2)

print(f"Conteggi sottesi alla gaussiana: {area_net:.3f} ±{sigma_area_net:.3f}")


# === Plot ===
c = ROOT.TCanvas("c", "Fit Gauss + Fondo", 800, 600)

# Imposta range X del disegno
H.GetXaxis().SetRangeUser(fit_min-100, fit_max+100)

H.Draw("E")  # Disegna istogramma con barre d'errore
fit_func.SetLineColor(ROOT.kRed)
fit_func.SetLineWidth(2)
fit_func.Draw("SAME")  # Sovrappone il fit

c.Update()
c.SaveAs("fit_gauss_parabola.png")


# === Calcolo del rateo di conteggio ===
time = 7 * 43000

rate_meas = area_net / time
print(f"Conteggi al secondo: {rate_meas:.3f}")

# === Calcolo del coefficiente di attenuazione ===
# --- Dati sperimentali ---
A_src = 188900         # Bq (fotoni/s)
err_A_src = 11647      # Bq (errore)

x = 1.0               # cm, spessore del bersaglio
BR = 0.903
fuga = 2.54 
diam_gate = 1.0 * 2.54        # cm

d_gate  = 16.0 + fuga         # cm
err_d_gate = 0.5              # cm (errore)

eff_spec = 0.33
err_eff_spec = 0.009780

eff_gate = 0.147
err_eff_gate = 0.01672  

# --- Angolo solido
beta = np.arctan((diam_gate/2) / d_gate)
Omega =  (1 - np.cos(beta))/2  #porzione di sfera

# PONENDO LO SPETTROMETRO A 25 CM IL CONO DELLA SORGENTE è CONTENUTO NELLA ZONA DELLA GEOMETRIA DEL DETECTOR
denominatore = A_src  * Omega * eff_gate * time * 2 * BR * eff_spec
lambda_val = - (1.0 / x) * np.log(area_net / denominatore)

print(f"Coefficiente di attenuazione λ = {lambda_val:.4f} cm⁻¹")
