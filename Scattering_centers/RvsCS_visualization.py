import matplotlib.pyplot as plt
import numpy as np
from iminuit import Minuit
from iminuit.cost import LeastSquares

file_path = "/mnt/c/Users/User/Desktop/info/Compton/Scattering_centers/"

def compton_energy(tetha):
    return 511 / (2 - np.cos(tetha))

def klein_nishina(theta):
    r = compton_energy(theta) / 511
    # The Klein-Nishina formula usually includes a factor of (r_e^2 / 2)
    # Here we only compute the angular dependent part, the constant factor will be fitted
    return (r ** 2) * (r + 1/r - np.sin(theta) ** 2)

thetas = np.linspace(np.radians(0), np.radians(120), 200) # Avoid 0 and pi for plotting stability

def errors(rate_errors, counts_errors, times):
    errors = []
    for i in range(len(rate_errors)):
        # Calculate the error for each rate
        rate_error = rate_errors[i]
        counts_error = counts_errors[i]
        time = times[i]

        error = np.sqrt((counts_error / time) ** 2 + (rate_error) ** 2)
        errors.append(error)
    return errors

angles = np.array([40, 50, 60, 70, 90, 110]) * np.pi / 180
measured_rates = np.array([0.03600, 0.0448, 0.03729, 0.0280, 0.03001, 0.01740])
measured_fit_errors = np.array([0.00031, 0.0004, 0.00035, 0.0002, 0.00022, 0.00017])
measured_counts_errors = np.sqrt(np.array([13692, 13371, 11165, 19239, 17892, 10994]))
measured_times = np.array([43000 * 8 + 36356, 43000 * 6 + 40184, 43000 * 6 + 31193 + 10254, 16 * 43000, 43000 * 13 + 37225, 43000 * 14 + 29739])
measured_errors = errors(measured_fit_errors, measured_counts_errors, measured_times)

# epsilon_gate = 0.476
simulated_rates = np.array([0.0871, 0.0619, 0.0456, 0.0365, 0.0256, 0.0242])
simulated_fit_errors = np.array([0.0017, 0.0010, 0.0011, 0.0007, 0.0005, 0.0005])
simulated_counts_errors = np.sqrt(np.array([34707, 37139, 17865, 21136, 22245, 13825]))
simulated_times = np.array([395946.67, 608513.63, 387603.03, 584773.04, 856244.45, 572827.43])
simulated_errors = errors(simulated_fit_errors, simulated_counts_errors, simulated_times)


# epsilon_gate = 0.166
# simulated_rates = np.array([0.0319, 0.0212, 0.0160, 0.0126, 0.0090, 0.0084])
# simulated_errors = np.array([0.0006, 0.0004, 0.0004, 0.0002, 0.0002, 0.0002])

# --- Fit Simulated Data ---
# Define the model function for the fit
def fit_model(theta, const):
    return const * klein_nishina(theta)

# Create the least squares cost function for simulated data
least_squares_sim = LeastSquares(angles, simulated_rates, simulated_errors, fit_model)

# Instantiate Minuit with an initial guess for the constant
# Initial guess can be estimated by Rate / KN(angle) for one point
initial_guess = simulated_rates[0] / klein_nishina(angles[0])
m_sim = Minuit(least_squares_sim, const=initial_guess)

# Run the minimization
m_sim.migrad()
m_sim.hesse() # Calculate uncertainties

# Get the fitted constant and its error
fitted_const = m_sim.values["const"]
fitted_const_err = m_sim.errors["const"]

print(f"Fitted constant for simulated data: {fitted_const:.5f} +/- {fitted_const_err:.5f}")
print(f"Fit Chi2/ndof: {m_sim.fval:.2f} / {m_sim.ndof}")

def constant(theta):
    r_gate = 1.27  #cm
    d_gate_source = 16 #cm
    alpha = np.arctan(r_gate / d_gate_source)
    Omega = 2 * np.pi * (1 - np.cos(alpha)) #sr

    r_sorg_targ = 5 #cm
    R = r_sorg_targ * np.tan(alpha)

    def A(theta, R):
        return np.pi * R ** 2 * np.cos(theta/2) 

    def V(theta, R, width):
        return A(theta, R) * width

    def N_c(theta, R, width):
        #For a Cu target
        rho = 8.96 #g/cm^3
        N_a = 6.022e23 #mol^-1
        Z = 29 #number of electrons
        MM = 63.546 #g/mol

        #Calculus of the number of atoms in the target
        N_c = rho * (N_a * Z / MM) * V(theta, R, width)
        return N_c

    def flux():
        S = 175000 #Bq
        BR = 0.903 * 2

        epsilon_gate = 0.16676 

        r_gate = 1.27  #cm
        A_gate = np.pi * r_gate ** 2 #cm^2

        # Calculus of the flux
        flux = (S * BR / A_gate) * (Omega / (4 * np.pi)) * epsilon_gate #cm^-2 s^-1
        return flux

    def epsilon_spettrometer(theta): 
        a = 0.19883
        b = -0.64438
        c = 0.11968

        return a * np.exp(-b * theta) + c

    r_e = 2.81794e-13 #cm

    return (r_e ** 2 / 2) * N_c(theta, R, 1) * flux() * epsilon_spettrometer(theta)/6

# --- Calculate Residuals using the fitted constant ---
residual_simulated = (simulated_rates - fitted_const * klein_nishina(angles))
residual_measured = (measured_rates - fitted_const * klein_nishina(angles))
residual_constant = (constant(thetas) * klein_nishina(thetas) - fitted_const * klein_nishina(thetas))

# --- Plotting ---
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True, gridspec_kw={'height_ratios': [3, 2]})

# Top plot: Data and Fit
ax1.errorbar(np.degrees(angles), simulated_rates, yerr=simulated_errors, fmt='o', label='Simulated', capsize=3, markersize=5)
ax1.errorbar(np.degrees(angles), measured_rates, yerr=measured_errors, fmt='o', label='Measured', capsize=3, markersize=5)
ax1.plot(np.degrees(thetas), fitted_const * klein_nishina(thetas),
         label=f'Fit (Simulated): c * KN\nc = {fitted_const:.3f} ± {fitted_const_err:.3f}',
         color='red')
ax1.plot(np.degrees(thetas), constant(thetas) * klein_nishina(thetas),
         label=r'$c(\theta) * KN$',
         color='blue')
ax1.set_ylabel('Rate (counts/s)')
ax1.set_title('Compton Scattering Rate vs Angle')
ax1.legend()
ax1.grid(True, linestyle=':')

# Bottom plot: Residuals
ax2.errorbar(np.degrees(angles), residual_simulated, fmt='o', label='Simulated Residual', capsize=3, markersize=5) # Error is 1 for standardized residuals
ax2.errorbar(np.degrees(angles), residual_measured, fmt='o', label='Measured Residual', capsize=3, markersize=5) # Error is 1 for standardized residuals
ax2.plot(np.degrees(thetas), residual_constant, label='Constant Residual', color='blue')
ax2.axhline(0, color='black', linestyle='--')
ax2.set_xlabel('Angle (degrees)')
ax2.set_ylabel('(Data - Fit) / Error')
ax2.legend()
ax2.grid(True, linestyle=':')

plt.tight_layout()
plt.savefig(file_path + 'compton_scattering_fit.png')
