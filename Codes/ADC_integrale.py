import ROOT
import MoraPyRoot as mpr
import matplotlib.pyplot as plt

file_path = "/mnt/c/Users/User/Desktop/info/Compton/ADC/"
# formatura = "preamplificata"
# formatura = "gaussiana-shaping"

# Funzione per leggere i dati da un file
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

def fill_TGraphError(formatura):
    amp_values, amp_errors = read_data(file_path + "data/formatura_" + formatura + "_amp(mV).txt")
    channel_values, channel_errors = read_data(file_path + "data/formatura_" + formatura + "_canale.txt")

    # Convert data in array for ROOT 
    import array
    amp_values_array = array.array('d', amp_values)
    channel_values_array = array.array('d', channel_values)
    amp_errors_array = array.array('d', amp_errors)
    channel_errors_array = array.array('d', channel_errors)

    graph_amp = ROOT.TGraphErrors(len(amp_values), channel_values_array, amp_values_array, channel_errors_array, amp_errors_array)

    return graph_amp

def extract_points_from_tgraph(graph):    
    points = []
    n_points = graph.GetN()  # Numero di punti nel TGraphErrors
    
    for i in range(n_points):
        x = graph.GetPointX(i)
        y = graph.GetPointY(i)
        err_x = graph.GetErrorX(i)
        err_y = graph.GetErrorY(i)
        points.append((float(x), float(y), float(err_x), float(err_y)))
    
    return points

def residuals(graph, fit_function, formatura):
    points = extract_points_from_tgraph(graph)
    
    res = []
    d = []

    for x, y, x_err, y_err in points: 
        res.append((y - fit_function.Eval(x)) / y_err)
        d.append(x)
            
    plt.figure(figsize=(10, 6))
    plt.errorbar(d, res, fmt='o', color='black')

    plt.axhline(0, color='red', linestyle='--', linewidth=1.5)

    plt.title("Residui formatura" + formatura)
    plt.xlabel(graph.GetXaxis().GetTitle())
    plt.ylabel(r'$\sigma$')
    # plt.ylim(-2, 2)
    plt.grid(True)
    plt.savefig(file_path + "formatura_" + formatura + "_residui.png")

coo1 = [0.1, 0.6, 0.45, 0.9]  
str1 = ["q", "m"] 

# Formatura gaussiana
formatura = "gaussiana"
T_gaus = fill_TGraphError(formatura)
f_gaus = ROOT.TF1("ADC", "pol1(0)", 0, 2000)
f_gaus.SetParameters(0, 0)
f_gaus.SetParameters(1, 0.2)

extreme = [0, 300, 0, 2000]
mpr.stampa_graph_fit(T_gaus, f_gaus, file_path +  "plots/ADC_integrale_formatura_" + formatura + ".png", "ADC integrale formatura " + formatura, "Canale", "Ampiezza (mV)", "AP", 300, 2000, 2, coo1, str1)
residuals(T_gaus, f_gaus, formatura)
mpr.stampa_graph_fit_range(T_gaus, f_gaus, extreme, file_path +  "plots/ADC_integrale_formatura_" + formatura + "_zoom.png", "ADC integrale formatura " + formatura, "Canale", "Ampiezza (mV)", "AP", 5, 300, 2, coo1, str1)

# Formatura preamplificataa
formatura = "preamplificata"
T_preamp = fill_TGraphError(formatura)
f_preamp = ROOT.TF1("ADC", "pol1(0)", 0, 2000)
f_preamp.SetParameters(0, 0)
f_preamp.SetParameters(1, 0.2)

extreme = [1600, 1900, 8000, 10000]
mpr.stampa_graph_fit(T_preamp, f_preamp, file_path +  "plots/ADC_integrale_formatura_" + formatura + ".png", "ADC integrale formatura " + formatura, "Canale", "Ampiezza (mV)", "AP", 950, 1600, 2, coo1, str1)
residuals(T_preamp, f_preamp, formatura)
mpr.stampa_graph_fit_range(T_preamp, f_preamp, extreme, file_path +  "plots/ADC_integrale_formatura_" + formatura + "_zoom.png", "ADC integrale formatura " + formatura, "Canale", "Ampiezza (mV)", "AP", 1600, 1900, 2, coo1, str1)

