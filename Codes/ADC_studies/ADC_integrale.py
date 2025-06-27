import ROOT
import matplotlib.pyplot as plt
import array

file_path = "/mnt/c/Users/ASUS/Desktop/WSL_Shared/Compton/Charatterization_elettronic_chain/ADC/"
formatura_list = ["gaussiana", "preamplificata", "gaussiana-shaping"]
title_list = ["Gaussian", "Preamplified", "Gaussian + shaping"]

# --- Funzione per leggere dati ---
def read_data(fileName):
    x_values, x_errors = [], []
    with open(fileName, 'r') as file:
        for line in file:
            parts = line.strip().split()
            if len(parts) == 2:
                x, dx = map(float, parts)
                x_values.append(x)
                x_errors.append(dx)
    return x_values, x_errors

# --- Costruisci TGraphErrors da file ---
def fill_TGraphError(formatura):
    amp_values, amp_errors = read_data(file_path + f"data/formatura_{formatura}_amp(mV).txt")
    channel_values, channel_errors = read_data(file_path + f"data/formatura_{formatura}_canale.txt")

    return ROOT.TGraphErrors(
        len(amp_values),
        array.array('d', channel_values),
        array.array('d', amp_values),
        array.array('d', channel_errors),
        array.array('d', amp_errors)
    )

# --- Estrai i punti da un TGraphErrors ---
def extract_points_from_tgraph(graph):    
    points = []
    for i in range(graph.GetN()):
        x = graph.GetPointX(i)
        y = graph.GetPointY(i)
        err_x = graph.GetErrorX(i)
        err_y = graph.GetErrorY(i)
        points.append((float(x), float(y), float(err_x), float(err_y)))
    return points

# --- Calcola e salva i residui ---
def plot_residuals(graph, fit_function, formatura):
    points = extract_points_from_tgraph(graph)
    x_vals, res_vals = [], []

    for x, y, _, y_err in points: 
        if y_err != 0:
            residual = (y - fit_function.Eval(x)) / y_err
        else:
            residual = 0
        x_vals.append(x)
        res_vals.append(residual)

    plt.figure(figsize=(10, 4))
    plt.axhline(0, color='red', linestyle='--', linewidth=1)
    plt.errorbar(x_vals, res_vals, fmt='o', color='black')
    plt.xlabel("Channel")
    plt.ylabel(r"Residual ($\sigma$)")
    plt.title(f"Residuals - {formatura}")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(file_path + f"plots/residui_{formatura}.png")
    plt.close()

# --- Inizializza grafico matplotlib unico ---
plt.figure(figsize=(16, 6))

# --- Ciclo sulle formature ---
for idx, formatura in enumerate(formatura_list):
    graph = fill_TGraphError(formatura)

    # Fit lineare
    fit_func = ROOT.TF1("fit", "pol1(0)", 0, 2000)
    fit_func.SetParameters(0, 0)
    fit_func.SetParameters(1, 0.2)
    graph.Fit(fit_func, "Q")  # 'Q' = quiet mode

    # Estrazione punti + plot
    points = extract_points_from_tgraph(graph)
    x_vals = [p[0] for p in points]
    y_vals = [p[1] for p in points]
    x_errs = [p[2] for p in points]
    y_errs = [p[3] for p in points]

    error = plt.errorbar(x_vals, y_vals, xerr=x_errs, yerr=y_errs, fmt='o', capsize=1)
    color = error.lines[0].get_color()
    # Plot funzione di fit
    x_fit = sorted(x_vals)
    y_fit = [fit_func.Eval(x) for x in x_fit]
    # Parametri del fit
    q = fit_func.GetParameter(0)
    m = fit_func.GetParameter(1)

    # Etichetta del fit con formula
    fit_label = f"Fit {title_list[idx]}: y = {m:.3f}x + {q:.1f}"
    plt.plot(x_fit, y_fit, '--', color=color, label=fit_label)


    # Salva residui per ogni formatura
    plot_residuals(graph, fit_func, formatura)

# --- Finalizza grafico matplotlib unico ---
plt.xlabel("Channel", fontsize=14)
plt.ylabel("Amplitude (mV)", fontsize=14)

plt.tick_params(labelsize=12)  # ingrandisce numeri sugli assi

plt.grid(True)
plt.xlim(0, 2000)
plt.legend(fontsize=16)
plt.tight_layout()
plt.savefig(file_path + "plots/ADC_integrale_TOT.png")

