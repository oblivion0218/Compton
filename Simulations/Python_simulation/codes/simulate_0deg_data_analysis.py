

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')  # allow saving figures without a display
import matplotlib.pyplot as plt
from tqdm import tqdm
from iminuit import Minuit
from iminuit.cost import LeastSquares


def gaussian(x, amp, mean, sigma):
    """Gaussian function."""
    return amp * np.exp(-0.5 * ((x - mean) / sigma) ** 2)


def linear(x, slope, intercept):
    """Linear background function."""
    return slope * x + intercept


def total_model(x, amp, mean, sigma, slope, intercept):
    """Combined Gaussian plus linear background."""
    return gaussian(x, amp, mean, sigma) + linear(x, slope, intercept)


def initial_guesses(x, y):
    """Estimate initial fit parameters from data."""
    idx_sorted = np.argsort(y)
    n_low = max(int(0.1 * len(y)), 1)
    intercept = np.median(y[idx_sorted[:n_low]])
    slope = 0.0
    amp = y.max() - intercept
    mean = x[y.argmax()]
    sigma = (x.max() - x.min()) / 6.0
    return amp, mean, sigma, slope, intercept


def compute_integrals(params, x_min, x_max):
    """Compute Gaussian and background integrals over [x_min, x_max]."""
    amp, mean, sigma, slope, intercept = params
    int_gauss = amp * sigma * np.sqrt(2 * np.pi)
    int_bg = slope / 2 * (x_max**2 - x_min**2) + intercept * (x_max - x_min)
    net = int_gauss - int_bg
    return int_gauss, int_bg, net


def fit_and_plot(file_path, plot_path, integrals_fh, n_bins=100):
    """Load energy list (skip header), histogram, fit model, save plot, and log net integral."""
    energies = np.loadtxt(file_path, skiprows=15)
    counts, edges = np.histogram(energies, bins=n_bins)
    x = (edges[:-1] + edges[1:]) / 2
    y = counts

    init = initial_guesses(x, y)

    # Define errors for least squares (sqrt of counts, min 1)
    errors = np.sqrt(y)
    errors[errors == 0] = 1.0

    # Set up least-squares cost: requires x, y, yerr, and model
    cost = LeastSquares(x, y, errors, total_model)

    # Initialize and configure Minuit
    m = Minuit(cost,
               amp=init[0], mean=init[1], sigma=init[2],
               slope=init[3], intercept=init[4])
    m.limits['sigma'] = (0, None)
    m.errordef = Minuit.LEAST_SQUARES
    m.migrad()

    params = (m.values['amp'], m.values['mean'], m.values['sigma'],
              m.values['slope'], m.values['intercept'])

    x_min, x_max = x.min(), x.max()
    _, _, net_int = compute_integrals(params, x_min, x_max)

    filename = os.path.basename(file_path)
    integrals_fh.write(f"{filename}\t{net_int:.6f}\n")

    # Plot
    y_fit = total_model(x, *params)
    y_gauss = gaussian(x, params[0], params[1], params[2])
    y_bg = linear(x, params[3], params[4])

    plt.figure(figsize=(8, 6))
    plt.bar(x, y, width=edges[1]-edges[0], alpha=0.6, label='Data')
    plt.plot(x, y_gauss, 'r--', label='Gaussian')
    plt.plot(x, y_bg, 'g--', label='Background')
    plt.plot(x, y_fit, 'k-', label='Total fit')
    plt.xlabel('Energy')
    plt.ylabel('Counts')
    plt.title(f"Fit for {filename}")
    plt.legend()
    plt.tight_layout()
    plt.savefig(plot_path)
    plt.close()


def main():
    data_dir = '/home/leonardo/Compton/Simulations/Python_simulation/simulated_events_NoEff/0_deg/data'
    plots_dir = '/home/leonardo/Compton/Simulations/Python_simulation/simulated_events_NoEff/0_deg/plots'

    os.makedirs(plots_dir, exist_ok=True)

    integrals_path = os.path.join(plots_dir, 'integrals.txt')
    with open(integrals_path, 'w') as fh:
        fh.write("# Filename\tNet Gaussian Integral (bg subtracted)\n")
        for fname in tqdm(sorted(os.listdir(data_dir))):
            if not fname.endswith('_detected_energies.txt'):
                continue
            file_path = os.path.join(data_dir, fname)
            plot_name = os.path.splitext(fname)[0] + '.png'
            plot_path = os.path.join(plots_dir, plot_name)
            try:
                fit_and_plot(file_path, plot_path, fh)
            except Exception as e:
                print(f"Error processing {fname}: {e}")

if __name__ == '__main__':
    main()