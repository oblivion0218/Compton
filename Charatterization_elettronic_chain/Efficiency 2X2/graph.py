import ROOT
import array
import numpy as np

# --------------------------------------------------------------------
# 1) DEFINIZIONE DEI DATI
# --------------------------------------------------------------------
# Abbiamo 3 set di dati: (x=Energia MeV, y=Efficienza)
#  - rivelatore da 2"
#  - rivelatore da 1.5"
#  - rivelatore da 3"

# --- Rivelatore 2" (7 punti) ---
x_2pollici = [
    0.31849513821826525,
    0.34839758515200353,
    0.5063605228428064,
    0.6578534128290872,
    0.7639818013788603,
    0.8357095062407706,
    1.1700423339274262
]

y_2pollici = [
    0.5530232692950782,
    0.513935011339386,
    0.3434113151787839,
    0.2491928618391422,
    0.20746631638582969,
    0.17917490984689505,
    0.12765545897224145
]

# --- Rivelatore 3" (7 punti) ---
x_3pollici = [
    0.31849513821826525,
    0.35101270017822783,
    0.5101613270630634,
    0.6627913412713632,
    0.7697163425286677,
    0.8357095062407706,
    1.1788248152624714,
]
y_1_5pollici = [
    0.3402790414078734,
    0.32207632355222093,
    0.20937604699537318,
    0.12765545897224145,
    0.1033982502308196,
    0.09608995471851453,
    0.06599675545070971,

]

# --- Rivelatore 1.5" (7 punti) ---
x_1_5pollici = [
    0.31849513821826525,
    0.3458019532572924,
    0.5063605228428064,
    0.6529522729442091,
    0.7639818013788603,
    0.8294833027836905,
    1.1788248152624714
]
y_3pollici = [
    0.6229821545978553,
    0.5789491301073059,
    0.42010843619294497,
    0.3191386474306036,
    0.2885399811814428,
    0.2706139612445929,
    0.2036989525676702,
]

Incertezza =  0.025 #L'HO SCELTA IO!

# --------------------------------------------------------------------
# 2) FUNZIONE DI FIT E FUNZIONE DI SUPPORTO PER FIT + GRAFICO
# --------------------------------------------------------------------
# Vogliamo f(x) = [0]*x^(-[1]) * exp(-[2]*x).
# Creiamo una funzione di comodo per:
# - Costruire il TGraph
# - Fare il Fit
# - Restituire (A, B, C) e gli oggetti TGraph, TF1

def fit_dataset(xvals, yvals, color, marker_style, legend_label):
    """Crea un TGraph, fitta con la funzione A*x^-B*exp(-C*x), restituisce:
       (graph, fit_func, A, B, C)."""
    n = len(xvals)
    ax = array.array("d", xvals)
    ay = array.array("d", yvals)
    ex = array.array("d", [0.0] * n)
    ey = array.array("d", [Incertezza * y for y in yvals])  
    
    graph = ROOT.TGraphErrors(n, ax, ay, ex ,ey)
    graph.SetMarkerColor(color)
    graph.SetLineColor(color)
    graph.SetMarkerStyle(marker_style)
    graph.SetMarkerSize(1.2)
    
    # Range degli x per la funzione di fit
    xmin, xmax = min(xvals), max(xvals)
    fit_func = ROOT.TF1("fit_"+legend_label, "[0]*pow(x, -[1])*exp(-[2]*x)+[3]", xmin, xmax)
    # Parametri iniziali (non troppo grandi)
    fit_func.SetParameters(1, 1, 1, 1)
    # Esegui il fit
    graph.Fit(fit_func, "QRS")  # Q=quiet, R=range, S=store
    
    A = fit_func.GetParameter(0)
    B = fit_func.GetParameter(1)
    C = fit_func.GetParameter(2)
    D = fit_func.GetParameter(3)
    
    return graph, fit_func, A, B, C, D


# --------------------------------------------------------------------
# 3) FITTIAMO I 3 SET DI DATI
# --------------------------------------------------------------------
# Assegniamo a ciascuno un colore e uno stile differente
graph_2, fit_2, A2, B2, C2 , D2 = fit_dataset(
    x_2pollici, y_2pollici,
    color=ROOT.kBlue, 
    marker_style=20,
    legend_label="2'' "
)

graph_1_5, fit_1_5, A1_5, B1_5, C1_5, D1_5 = fit_dataset(
    x_1_5pollici, y_1_5pollici,
    color=ROOT.kRed,
    marker_style=21,
    legend_label="1_5''"
)

graph_3, fit_3, A3, B3, C3, D3 = fit_dataset(
    x_3pollici, y_3pollici,
    color=ROOT.kGreen+2,
    marker_style=22,
    legend_label="3''"
)


# --------------------------------------------------------------------
# 4) INTERPOLAZIONE/ESTRAPOLAZIONE DEI PARAMETRI PER 1" 
# --------------------------------------------------------------------
# Abbiamo (dim in pollici) -> (A, B, C):
# Facciamo un fit lineare di A(d), B(d), C(d) in funzione di d (d=1.5,2,3)
# e valutiamo a d=1 pollice.

fit_511_params = {}

def interpolate_parameter(d, Y, i):
 # Costruzione di TGraphErrors con errore 5% su Y
    n = len(d)
    ax = array.array("d", d)
    ay = array.array("d", Y)
    ex = array.array("d", [0.0] * n)
    ey = array.array("d", [Incertezza * y for y in Y])  # 5% errore relativo

    g = ROOT.TGraphErrors(n, ax, ay, ex, ey)

    # Fit function: f(d) = p3 + p0 * exp(p1 * d - p2)
    f_exp = ROOT.TF1("f_exp", " [0]*x + [1]", 0, 5)
    f_exp.SetParameters(1, 1)  # parametri iniziali

    g.Fit(f_exp, "Q")  # Quiet fit

    # Estrai parametri e incertezze
    p = [f_exp.GetParameter(k) for k in range(2)]
    dp = [f_exp.GetParError(k) for k in range(2)]

    # Salva i parametri se i == 2 (511 keV)
    if i == 2:
        fit_511_params["params"] = p
        fit_511_params["errors"] = dp
        print("Parametri e errori salvati per 511 keV (1 pollice).")

    # Plot
    c_fit = ROOT.TCanvas("c_fit", f"Fit Param {i}", 800, 600)
    g.SetTitle("Fit dei parametri in funzione della dimensione; Dimensione (pollici); Valore del parametro")
    g.SetMarkerStyle(20)
    g.SetMarkerColor(ROOT.kBlue)
    g.Draw("AP")
    f_exp.SetLineColor(ROOT.kRed)
    f_exp.Draw("same")
    c_fit.SaveAs(f"fit_parameter{i}.png")

    # Valore interpolato a d = 1
    return  p[0] * 1 + p[1]  

diams = [1.5, 2.0, 3.0]
y_1 = []

for i in range(7):
    y_1.append(interpolate_parameter(diams, [y_1_5pollici[i], y_2pollici[i], y_3pollici[i]], i))

graph_1, fit_1, A1, B1, C1, D1 = fit_dataset(
    x_2pollici, y_1,
    color=ROOT.kMagenta, 
    marker_style=24,
    legend_label="1'' "
)

def predict_efficiency_511 (fit_params, d=1):
    # Estrai parametri e incertezze
    p0, p1 = fit_params['params']
    dp0, dp1 = fit_params['errors']

    print( "P0:", p0, "±", dp0)
    print( "P1:", p1, "±", dp1)
    

    eff = p0 * d + p1
    delta_eff = np.sqrt((dp0 * d)**2 + (dp1)**2)

    return eff, delta_eff

eff_511, err_511 = predict_efficiency_511(fit_511_params)
print(f"\n Efficienza a 511 keV (1 pollice): {eff_511:.5f} ± {err_511:.5f} \n")
print("ATTENZIONE L'ERRORE SOPRA è DATO DA GLI ERRORI DEI PARAMETRI DEL FIT, NON DA UN ERRORE STATISTICO!")
print("Se impongo come nei casi precedenti un errore del ", Incertezza*100 , "% su ogni punto, l'errore finale sarà: " , eff_511 *  Incertezza)

# --------------------------------------------------------------------
# 5) DISEGNO SU UN UNICO CANVAS + LEGENDA
# --------------------------------------------------------------------
c = ROOT.TCanvas("c", "Efficienza vs Energia", 1000, 700)
c.SetGrid()
#c.SetLogx()
#c.SetLogy()

# Per disegnare più TGraph in un unico asse, usiamo un TMultiGraph
mg = ROOT.TMultiGraph()

mg.Add(graph_2, "P")       # "P" = draw markers
mg.Add(graph_1_5, "P")
mg.Add(graph_3, "P")
mg.Add(graph_1,"P")

# Disegniamo il MultiGraph
mg.SetTitle("Confronto rivelatori e fit; Gamma energy (MeV); Intrinsic peak efficiency")
#mg.GetXaxis().SetLimits(0, 3)  # Imposta il range dell'asse X
mg.Draw("A")

# Ora disegniamo le curve di fit di ciascuno con "same"
fit_2.SetLineColor(ROOT.kBlue)
fit_1_5.SetLineColor(ROOT.kRed)
fit_3.SetLineColor(ROOT.kGreen+2)
fit_1.SetLineColor(ROOT.kMagenta)
fit_1.SetLineStyle(2)

fit_2.Draw("same")
fit_1_5.Draw("same")
fit_3.Draw("same")
fit_1.Draw("same")


# Creiamo la leggenda
legend = ROOT.TLegend(0.65, 0.65, 0.88, 0.88)
legend.SetBorderSize(0)
legend.SetFillStyle(0)

# Aggiungiamo le entry
legend.AddEntry(graph_1_5,f"1.5 inches ", "lp")
legend.AddEntry(graph_2,  f"2 inchesi", "lp")
legend.AddEntry(graph_3,  f"3 inches", "lp")
legend.AddEntry(graph_1, f"Prev. 1 inch", "lp")

legend.SetTextSize(0.03)
legend.Draw()

c.Update()
c.SaveAs("Efficienze.png")


#-------------- STAMPA VALORI A SCHERMO --------

print("\nParametri del fit per il rivelatore da 2 pollici (A*x^-B*exp(-C*x)+D):")
print(f"  A = {fit_2.GetParameter(0):.5f} ± {fit_2.GetParError(0):.5f}")
print(f"  B = {fit_2.GetParameter(1):.5f} ± {fit_2.GetParError(1):.5f}")
print(f"  C = {fit_2.GetParameter(2):.5f} ± {fit_2.GetParError(2):.5f}")
print(f"  D = {fit_2.GetParameter(3):.5f} ± {fit_2.GetParError(3):.5f}\n")

print("Previsioni efficienza per 2'' usando il fit:")
print(f"x = 0.511 MeV → efficienza = {fit_2.Eval(0.511):.5f} ")
print(f"x = 1.173 MeV → efficienza = {fit_2.Eval(1.173):.5f} ")
print(f"x = 1.274 MeV → efficienza = {fit_2.Eval(1.274):.5f} ")
print(f"x = 1.332 MeV → efficienza = {fit_2.Eval(1.332):.5f} \n")

ang_test = [15, 35, 40, 50, 60, 70, 90, 110]  # Angoli in gradi
ang = [ang * 0.01745 for ang in ang_test]  # Conversione in radianti

en_compton = []

# Calcolare l'energia di Compton per ogni angolo
for i in range(len(ang)):
    en_compton.append(511 / (2 - np.cos(ang[i])))

print("Previsioni efficienza per 2'' usando il fit:")
for i, x in enumerate(en_compton):  
    y_pred = fit_2.Eval(x/1000)  
    print(f"Angolo = {ang_test[i]}° ----- x = {x:.3f} KeV → efficienza = {y_pred:.5f}")