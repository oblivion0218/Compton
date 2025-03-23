import numpy as np
import matplotlib.pyplot as plt
from iminuit import Minuit
from scipy.stats import norm

# === Leggi i dati ===
def read_data(filename):
    data = []
    with open(filename, 'r') as f:
        for line in f:
            if line.strip() and not line.startswith("Canale"):
                try:
                    _, _, eventi = map(int, line.split())
                    data.append(eventi)
                except:
                    continue
    return np.array(data)

filename = "tot_15.txt"
data = read_data(filename)
x_vals = np.arange(len(data))
mask = (x_vals >= 1100) & (x_vals <= 1800)
x = x_vals[mask]
y = data[mask]
y_err = np.sqrt(np.maximum(y, 1))

# === Funzione totale (doppia gaussiana + parabola) ===
def model(x, A1, mu1, sigma1, A2, mu2, sigma2, a0, a1, a2):
    g1 = A1 * np.exp(-0.5 * ((x - mu1) / sigma1)**2)
    g2 = A2 * np.exp(-0.5 * ((x - mu2) / sigma2)**2)
    bkg = a0 + a1 * x + a2 * x**2
    return g1 + g2 + bkg

# === Loop con fit ===
for i in range(20):
    mu2_fixed = 1490 + i

    # Definizione funzione chi2 con mu2 fisso
    def chi2(A1, mu1, sigma1, A2, sigma2, a0, a1, a2):
        y_fit = model(x, A1, mu1, sigma1, A2, mu2_fixed, sigma2, a0, a1, a2)
        return np.sum(((y - y_fit) / y_err) ** 2)

    # Fit con Minuit
    m = Minuit(chi2,
               A1=150, mu1=1450, sigma1=100,
               A2=200, sigma2=100,
               a0=0, a1=0, a2=0)

    m.limits["A1"] = (100, 200)
    m.limits["mu1"] = (1400, 1500)
    m.limits["sigma1"] = (10, 120)
    m.limits["A2"] = (50, 200)
    m.limits["sigma2"] = (10, 200)
    m.errordef = Minuit.LEAST_SQUARES
    m.migrad()

    # === Plot ===
    fig, ax = plt.subplots()
    ax.errorbar(x, y, yerr=y_err, fmt='.', label='Dati')

    x_fit = np.linspace(1100, 1800, 1000)
    y_fit = model(x_fit, *m.values.values(), mu2=mu2_fixed)
    ax.plot(x_fit, y_fit, 'r-', label='Fit Totale')

    ax.plot(x_fit, m.values["A1"] * np.exp(-0.5 * ((x_fit - m.values["mu1"]) / m.values["sigma1"])**2),
            '--', label='Compton 15°', color='green')

    ax.plot(x_fit, m.values["A2"] * np.exp(-0.5 * ((x_fit - mu2_fixed) / m.values["sigma2"])**2),
            '--', label='511 keV', color='blue')

    ax.plot(x_fit, m.values["a0"] + m.values["a1"] * x_fit + m.values["a2"] * x_fit**2,
            '--', label='Fondo', color='orange')

    ax.set_title(f"Fit Istogramma {mu2_fixed}")
    ax.set_xlabel("Canale")
    ax.set_ylabel("Conteggi")
    ax.legend()
    ax.grid(True)

    plt.tight_layout()
    plt.savefig(f"Fit_Istogramma_{mu2_fixed}.png")
    plt.close()
