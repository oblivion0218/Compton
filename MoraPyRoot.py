import ROOT
import numpy as np
import matplotlib.pyplot as plt

class ScientificNotation:
    def __init__(self, n, exp):
        self.n = n
        self.exp = exp

def exponential(n):
    i = 0
    if n < 1:
        while n < 1:
            n *= 10
            i -= 1
    elif n >= 10:
        while n >= 10:
            n /= 10
            i += 1
    return ScientificNotation(n, i)

# Serve per importare un oggetto da un root file
def import_Tobject(file_name, object_name):
    file = ROOT.TFile(file_name)
    if file.IsZombie():
        raise Exception("Error opening the ROOT file")
    
    obj = file.Get(object_name)
    if obj is None:
        raise Exception(f"Error finding object: {object_name}")
    
    return obj.Clone(object_name)

#-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
# Funzioni per fit
#-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

def fit(point, function, n_parameters, option, precision, min_val=None, max_val=None, cov_mat=False):
    fit_result = point.Fit(function, option, "", min_val, max_val) if min_val and max_val else point.Fit(function, option)
    
    print("\n\nFit result:", fit_result.IsValid())
    print("\n")
    for i in range(n_parameters):
        print(f"par {i}:\t{function.GetParameter(i)} ± {function.GetParError(i)}")
    
    print("\nDegrees of freedom:", fit_result.Ndf())
    print("Chi2:", fit_result.Chi2())
    print("Probability:", fit_result.Prob(), "\t", fit_result.Prob() * 100, "%\n\n")
    
    if cov_mat:
        fit_result.PrintCovMatrix(ROOT.std.cout)

def stampa_graph_fit(point, function, destination_png, graph_name, x_axis_name, y_axis_name, graphic_option, min_val = None, max_val = None, n_parameters = 0, pave_coordinates = None, pave_entries = None):
    canvas = ROOT.TCanvas()
    
    point.Draw(graphic_option)
    point.SetTitle(graph_name)
    point.GetXaxis().SetTitle(x_axis_name)
    point.GetYaxis().SetTitle(y_axis_name)

    fit_result = point.Fit(function, "S", "", min_val, max_val) if min_val and max_val else point.Fit(function, "S")

    if pave_coordinates != None and pave_entries != None : 
        text_box = ROOT.TPaveText(pave_coordinates[0], pave_coordinates[1], pave_coordinates[2], pave_coordinates[3], "NDC")

        chi2 = exponential(fit_result.Chi2())
        par_errors = [exponential(function.GetParError(i)) for i in range(n_parameters)]

        text_box.SetFillColor(0)
        text_box.SetTextAlign(12)
        text_box.AddText("Fit result:")

        if abs(chi2.exp) < 3:
            text_box.AddText(f"#chi^{{2}}/dof = {fit_result.Chi2():.3f}/{fit_result.Ndf()} = {fit_result.Chi2() / fit_result.Ndf():.3f}")
        else:
            text_box.AddText(f"#chi^{{2}}/dof = {chi2.n:.2f} * 10^{{{chi2.exp}}}/{fit_result.Ndf()} = {fit_result.Chi2() / fit_result.Ndf():.2f} * 10^{{{exponential(fit_result.Chi2() / fit_result.Ndf()).exp}}}")

        for i in range(n_parameters):
            if abs(par_errors[i].exp) < 3:
                text_box.AddText(f"{pave_entries[i]} = {function.GetParameter(i):.3f} #pm {function.GetParError(i):.3f}")
            else:
                text_box.AddText(f"{pave_entries[i]} = ({function.GetParameter(i) * 10 ** -par_errors[i].exp:.3f} #pm {par_errors[i].n:.3f}) * 10^{{{par_errors[i].exp}}}")

        text_box.Draw()

        canvas.Print(destination_png, "png")
        del text_box
    
    else:            
        canvas.Print(destination_png, "png")

# SetRange
def set_range(point, ax, coordinates):
    if ax == "x" or ax == "X":
        point.GetXaxis().SetRangeUser(coordinates[0], coordinates[1])
    elif ax == "y" or ax == "Y":
        point.GetYaxis().SetRangeUser(coordinates[0], coordinates[1])
    elif ax == "xy" or ax == "XY":
        point.GetXaxis().SetRangeUser(coordinates[0], coordinates[1])
        point.GetYaxis().SetRangeUser(coordinates[2], coordinates[3])
    elif ax == "yx" or ax == "YX":
        point.GetYaxis().SetRangeUser(coordinates[0], coordinates[1])
        point.GetXaxis().SetRangeUser(coordinates[2], coordinates[3])
    
# extreme_graph = ("xy", [x_min, x_max, y_muin, y_max])
def stampa_graph_fit_range(point, function, extreme_graph, destination_png, graph_name, x_axis_name, y_axis_name, graphic_option, min_val = None, max_val = None, n_parameters = 0, pave_coordinates = None, pave_entries = None):
    canvas = ROOT.TCanvas()
    
    point.Draw(graphic_option)
    point.SetTitle(graph_name)
    point.GetXaxis().SetTitle(x_axis_name)
    point.GetYaxis().SetTitle(y_axis_name)

    fit_result = point.Fit(function, "S", "", min_val, max_val) if min_val and max_val else point.Fit(function, "S")

    if pave_coordinates != None and pave_entries != None : 
        text_box = ROOT.TPaveText(pave_coordinates[0], pave_coordinates[1], pave_coordinates[2], pave_coordinates[3], "NDC")

        chi2 = exponential(fit_result.Chi2())
        par_errors = [exponential(function.GetParError(i)) for i in range(n_parameters)]

        text_box.SetFillColor(0)
        text_box.SetTextAlign(12)
        text_box.AddText("Fit result:")

        if abs(chi2.exp) < 3:
            text_box.AddText(f"#chi^{{2}}/dof = {fit_result.Chi2():.3f}/{fit_result.Ndf()} = {fit_result.Chi2() / fit_result.Ndf():.3f}")
        else:
            text_box.AddText(f"#chi^{{2}}/dof = {chi2.n:.2f} * 10^{{{chi2.exp}}}/{fit_result.Ndf()} = {fit_result.Chi2() / fit_result.Ndf():.2f} * 10^{{{exponential(fit_result.Chi2() / fit_result.Ndf()).exp}}}")

        for i in range(n_parameters):
            if abs(par_errors[i].exp) < 3:
                text_box.AddText(f"{pave_entries[i]} = {function.GetParameter(i):.3f} #pm {function.GetParError(i):.3f}")
            else:
                text_box.AddText(f"{pave_entries[i]} = ({function.GetParameter(i) * 10 ** -par_errors[i].exp:.3f} #pm {par_errors[i].n:.3f}) * 10^{{{par_errors[i].exp}}}")

        text_box.Draw()

        set_range(point, extreme_graph[0], extreme_graph[1])

        canvas.Print(destination_png, "png")
        del text_box
    
    else:            
        set_range(point, extreme_graph[0], extreme_graph[1])
                
        canvas.Print(destination_png, "png")

#-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
# Funzioni per stampare velocemente con matplotlib un oggetto root
#-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

def plot_hist_MPL(hist, fileNamePNG, grid = True):
    # Estrai il numero di bin, i contenuti e i bordi del TH1D
    n_bins = hist.GetNbinsX()
    bin_contents = [hist.GetBinContent(i+1) for i in range(n_bins)]
    bin_edges = [hist.GetBinLowEdge(i+1) for i in range(n_bins+1)]
    
    # Plotta usando matplotlib
    plt.figure(figsize=(8, 6))
    plt.hist(bin_edges[:-1], bins=bin_edges, weights=bin_contents, histtype='step')
    plt.xlabel(hist.GetXaxis().GetTitle())
    plt.ylabel(hist.GetYaxis().GetTitle())
    plt.grid(grid)
    plt.savefig(fileNamePNG)

def plot_TGraph_MPL(graph, fileNamePNG, grid=True):
    # Estrai i punti (x, y) dal TGraph
    n_points = graph.GetN()
    x_values = [graph.GetPointX(i) for i in range(n_points)]
    y_values = [graph.GetPointY(i) for i in range(n_points)]
    
    # Plotta con matplotlib
    plt.figure(figsize=(8, 6))
    plt.plot(x_values, y_values, 'o-')
    plt.xlabel(graph.GetXaxis().GetTitle())
    plt.ylabel(graph.GetYaxis().GetTitle())
    plt.grid(grid)
    plt.savefig(fileNamePNG)

def plot_TGraphErrors_MPL(graph, fileNamePNG, grid=True):
    # Estrai i punti (x, y) e le incertezze dal TGraphErrors
    n_points = graph.GetN()
    x_values = [graph.GetPointX(i) for i in range(n_points)]
    y_values = [graph.GetPointY(i) for i in range(n_points)]
    x_errors = [graph.GetErrorX(i) for i in range(n_points)]
    y_errors = [graph.GetErrorY(i) for i in range(n_points)]
    
    # Plotta con matplotlib
    plt.figure(figsize=(8, 6))
    plt.errorbar(x_values, y_values, xerr=x_errors, yerr=y_errors, fmt='o')
    plt.xlabel(graph.GetXaxis().GetTitle())
    plt.ylabel(graph.GetYaxis().GetTitle())
    plt.grid(grid)
    plt.savefig(fileNamePNG)

def plot_TF1_MPL(function, x_min, x_max, fileNamePNG, grid=True, n_points=1000):
    # Genera i punti x e calcola y
    x_values = np.linspace(x_min, x_max, n_points)
    y_values = [function.Eval(x) for x in x_values]
    
    # Plotta con matplotlib
    plt.figure(figsize=(8, 6))
    plt.plot(x_values, y_values, '-')
    plt.xlabel("x")
    plt.ylabel("f(x)")
    plt.grid(grid)
    plt.savefig(fileNamePNG)

