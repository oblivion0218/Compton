import ROOT
import numpy as np
import matplotlib.pyplot as plt

class ScientificNotation:
    """
    Represents a number in scientific notation.
    
    :param n: The coefficient of the number.
    :param exp: The exponent of the number.
    """
    def __init__(self, n, exp):
        self.n = n
        self.exp = exp

def exponential(n):
    """
    Converts a number to scientific notation.
    
    :param n: The number to convert.
    :return: A ScientificNotation object.
    """
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

def import_Tobject(file_name, object_name):
    """
    Imports an object from a ROOT file.
    
    :param file_name: Path to the ROOT file.
    :param object_name: Name of the object to import.
    :return: A cloned instance of the requested object.
    """
    file = ROOT.TFile(file_name)
    if file.IsZombie():
        raise Exception("Error opening the ROOT file")
    
    obj = file.Get(object_name)
    if obj is None:
        raise Exception(f"Error finding object: {object_name}")
    
    return obj.Clone(object_name)

#-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
# Function to fit a ROOT object with a given function
#-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

def fit(point, function, n_parameters, option, precision, min_val=None, max_val=None, cov_mat=False):
    """
    Performs a fit on a given dataset using a specified function.
    
    :param point: ROOT object containing the data to be fitted.
    :param function: ROOT function used for the fit.
    :param n_parameters: Number of parameters in the fit function.
    :param option: Fit options string for ROOT.
    :param precision: Fit precision.
    :param min_val: Minimum x-range for the fit (optional).
    :param max_val: Maximum x-range for the fit (optional).
    :param cov_mat: If True, prints the covariance matrix of the fit.
    """
    fit_result = point.Fit(function, option, "", min_val, max_val) if min_val and max_val else point.Fit(function, option)
    
    print("\n\nFit result:", fit_result.IsValid())
    for i in range(n_parameters):
        print(f"par {i}:\t{function.GetParameter(i)} ± {function.GetParError(i)}")
    
    print("\nDegrees of freedom:", fit_result.Ndf())
    print("Chi2:", fit_result.Chi2())
    print("Probability:", fit_result.Prob(), "\t", fit_result.Prob() * 100, "%\n\n")
    
    if cov_mat:
        fit_result.PrintCovMatrix(ROOT.std.cout)

def stampa_graph_fit(point, function, destination_png, graph_name, x_axis_name, y_axis_name, graphic_option, min_val = None, max_val = None, n_parameters = 0, pave_coordinates = None, pave_entries = None):
    """
    Plots a ROOT graph and fits it with a given function.

    :param point: ROOT object containing the data to be plotted.
    :param function: ROOT function used for the fit.
    :param destination_png: Path to save the plot.
    :param graph_name: Title of the plot.
    :param x_axis_name: Label for the x-axis.
    :param y_axis_name: Label for the y-axis.
    :param graphic_option: Options for the plot.
    :param min_val: Minimum x-value for the fit range (optional).
    :param max_val: Maximum x-value for the fit range (optional).
    :param n_parameters: Number of parameters in the fit function.
    :param pave_coordinates: Coordinates for the fit results box (optional).
    :param pave_entries: Labels for the fit results box (optional).
    """
    canvas = ROOT.TCanvas()
    
    point.Draw(graphic_option)
    point.SetTitle(graph_name)
    point.GetXaxis().SetTitle(x_axis_name)
    point.GetYaxis().SetTitle(y_axis_name)

    fit_result = point.Fit(function, "S", "", min_val, max_val) if min_val and max_val else point.Fit(function, "S")

    if pave_coordinates != None and pave_entries != None : 
        text_box = ROOT.TPaveText(pave_coordinates[0], pave_coordinates[1], pave_coordinates[2], pave_coordinates[3], "NDC")

        # chi2 = exponential(fit_result.Chi2())
        par_errors = [exponential(function.GetParError(i)) for i in range(n_parameters)]

        text_box.SetFillColor(0)
        text_box.SetTextAlign(12)
        text_box.AddText("Fit result:")

        # if abs(chi2.exp) < 3:
        #     text_box.AddText(f"#chi^{{2}}/dof = {fit_result.Chi2():.3f}/{fit_result.Ndf()} = {fit_result.Chi2() / fit_result.Ndf():.3f}")
        # else:
        #     text_box.AddText(f"#chi^{{2}}/dof = {chi2.n:.2f} * 10^{{{chi2.exp}}}/{fit_result.Ndf()} = {fit_result.Chi2() / fit_result.Ndf():.2f} * 10^{{{exponential(fit_result.Chi2() / fit_result.Ndf()).exp}}}")

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

    return fit_result

def set_range(point, ax, coordinates):
    """
    Sets the range of the x or y axis of a ROOT graph.
    
    :param point: ROOT object whose axis range is being modified.
    :param ax: Axis to modify ('x', 'y', 'xy', or 'yx').
    :param coordinates: List of min and max values for the axis.
    """
    if ax.lower() == "x":
        point.GetXaxis().SetRangeUser(coordinates[0], coordinates[1])
    elif ax.lower() == "y":
        point.GetYaxis().SetRangeUser(coordinates[0], coordinates[1])
    elif ax.lower() in ["xy", "yx"]:
        point.GetXaxis().SetRangeUser(coordinates[0], coordinates[1])
        point.GetYaxis().SetRangeUser(coordinates[2], coordinates[3])

# extreme_graph = ("xy", [x_min, x_max, y_muin, y_max])
def stampa_graph_fit_range(point, function, extreme_graph, destination_png, graph_name, x_axis_name, y_axis_name, graphic_option, min_val = None, max_val = None, n_parameters = 0, pave_coordinates = None, pave_entries = None):
    """
    Plots a ROOT graph and fits it with a given function, setting the axis range.

    :param point: ROOT object containing the data to be plotted.
    :param function: ROOT function used for the fit.
    :param extreme_graph: Tuple containing the axis to modify and the range values.
    :param destination_png: Path to save the plot.
    :param graph_name: Title of the plot.
    :param x_axis_name: Label for the x-axis.
    :param y_axis_name: Label for the y-axis.
    :param graphic_option: Options for the plot.
    :param min_val: Minimum x-value for the fit range (optional).
    :param max_val: Maximum x-value for the fit range (optional).
    :param n_parameters: Number of parameters in the fit function.
    :param pave_coordinates: Coordinates for the fit results box (optional).
    :param pave_entries: Labels for the fit results box (optional).   
    """    
    canvas = ROOT.TCanvas()
    
    point.Draw(graphic_option)
    point.SetTitle(graph_name)
    point.GetXaxis().SetTitle(x_axis_name)
    point.GetYaxis().SetTitle(y_axis_name)

    fit_result = point.Fit(function, "S", "", min_val, max_val) if min_val and max_val else point.Fit(function, "S")

    if pave_coordinates != None and pave_entries != None : 
        text_box = ROOT.TPaveText(pave_coordinates[0], pave_coordinates[1], pave_coordinates[2], pave_coordinates[3], "NDC")

        # chi2 = exponential(fit_result.Chi2())
        par_errors = [exponential(function.GetParError(i)) for i in range(n_parameters)]

        text_box.SetFillColor(0)
        text_box.SetTextAlign(12)
        text_box.AddText("Fit result:")

        # if abs(chi2.exp) < 3:
        #     text_box.AddText(f"#chi^{{2}}/dof = {fit_result.Chi2():.3f}/{fit_result.Ndf()} = {fit_result.Chi2() / fit_result.Ndf():.3f}")
        # else:
        #     text_box.AddText(f"#chi^{{2}}/dof = {chi2.n:.2f} * 10^{{{chi2.exp}}}/{fit_result.Ndf()} = {fit_result.Chi2() / fit_result.Ndf():.2f} * 10^{{{exponential(fit_result.Chi2() / fit_result.Ndf()).exp}}}")

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

    return fit_result


#-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
# Function to plot with matplotlib a ROOT object
#-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

def plot_hist_MPL(hist, fileNamePNG, grid=True):
    """
    Plots a ROOT histogram using Matplotlib.
    
    :param hist: ROOT histogram object.
    :param fileNamePNG: Path to save the plot.
    :param grid: Whether to display the grid (default: True).
    """
    n_bins = hist.GetNbinsX()
    bin_contents = [hist.GetBinContent(i+1) for i in range(n_bins)]
    bin_edges = [hist.GetBinLowEdge(i+1) for i in range(n_bins+1)]
    
    plt.figure(figsize=(8, 6))
    plt.hist(bin_edges[:-1], bins=bin_edges, weights=bin_contents, histtype='step')
    plt.xlabel(hist.GetXaxis().GetTitle())
    plt.ylabel(hist.GetYaxis().GetTitle())
    plt.grid(grid)
    plt.savefig(fileNamePNG)

def plot_TGraph_MPL(graph, fileNamePNG, grid=True):
    """
    Plots a ROOT TGraph using Matplotlib.
    
    :param graph: ROOT TGraph object.
    :param fileNamePNG: Path to save the plot.
    :param grid: Whether to display the grid (default: True).
    """
    n_points = graph.GetN()
    x_values = [graph.GetPointX(i) for i in range(n_points)]
    y_values = [graph.GetPointY(i) for i in range(n_points)]
    
    plt.figure(figsize=(8, 6))
    plt.plot(x_values, y_values, 'o-')
    plt.xlabel(graph.GetXaxis().GetTitle())
    plt.ylabel(graph.GetYaxis().GetTitle())
    plt.grid(grid)
    plt.savefig(fileNamePNG)

def plot_TGraphErrors_MPL(graph, fileNamePNG, grid=True):
    """
    Plots a ROOT TGraphErrors using Matplotlib.
    
    :param graph: ROOT TGraph object.
    :param fileNamePNG: Path to save the plot.
    :param grid: Whether to display the grid (default: True).
    """
    n_points = graph.GetN()
    x_values = [graph.GetPointX(i) for i in range(n_points)]
    y_values = [graph.GetPointY(i) for i in range(n_points)]
    x_errors = [graph.GetErrorX(i) for i in range(n_points)]
    y_errors = [graph.GetErrorY(i) for i in range(n_points)]
    
    plt.figure(figsize=(8, 6))
    plt.errorbar(x_values, y_values, xerr=x_errors, yerr=y_errors, fmt='o')
    plt.xlabel(graph.GetXaxis().GetTitle())
    plt.ylabel(graph.GetYaxis().GetTitle())
    plt.grid(grid)
    plt.savefig(fileNamePNG)

def plot_TF1_MPL(function, x_min, x_max, fileNamePNG, grid=True, n_points=1000):
    """
    Plots a ROOT TF1 function using Matplotlib.
    
    :param function: ROOT TF1 function to plot.
    :param x_min: Minimum x-value.
    :param x_max: Maximum x-value.
    :param fileNamePNG: Path to save the plot.
    :param grid: Whether to display the grid (default: True).
    :param n_points: Number of points to sample in the plot.
    """
    x_values = np.linspace(x_min, x_max, n_points)
    y_values = [function.Eval(x) for x in x_values]
    
    plt.figure(figsize=(8, 6))
    plt.plot(x_values, y_values, '-')
    plt.xlabel("x")
    plt.ylabel("f(x)")
    plt.grid(grid)
    plt.savefig(fileNamePNG)

