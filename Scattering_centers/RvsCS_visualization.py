import matplotlib.pyplot as plt
import numpy as np

def compton_energy(tetha):
    return 511 / (2 - np.cos(tetha))

def klein_nishina(theta):
    r = compton_energy(theta) / 511
    return r ** 2 * (r + 1/r - np.sin(theta) ** 2)

thetas = np.linspace(0, np.pi, 100)

angles = np.array([40, 60, 70, 90, 110]) * np.pi / 180
measured_rates = np.array([0.03600, 0.03729, 0.0515, 0.03001, 0.01740])
measured_errors = np.array([0.00031, 0.00035, 0.0004, 0.00022, 0.00017])
simulated_rates = np.array([0.0871, 0.0456, 0.0365, 0.0256, 0.0242])
simulated_errors = np.array([0.0017, 0.0011, 0.0007, 0.0005, 0.0005])

# const = 1.7737663904882119e-43
const = 0.14327846609027658 * 0.5
residual_simulated = (simulated_rates - const * klein_nishina(angles)) / (const * klein_nishina(angles))
residual_measured = (measured_rates - const * klein_nishina(angles)) / (const * klein_nishina(angles))

# Plotting upper the klein-nishina and lower the residual
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
ax1.errorbar(angles, simulated_rates, yerr=simulated_errors, fmt='o', label='Simulated')
ax1.errorbar(angles, measured_rates, yerr=measured_errors, fmt='o', label='Measured')
ax1.plot(thetas, const * klein_nishina(thetas), label='Klein-Nishina', color='red')
ax1.set_ylabel('Rate (counts/s)')
ax1.legend()
ax2.errorbar(angles, residual_simulated, yerr=simulated_errors / simulated_rates, fmt='o', label='Simulated Residual')
ax2.errorbar(angles, residual_measured, yerr=measured_errors / measured_rates, fmt='o', label='Measured Residual')
ax2.axhline(0, color='black', linestyle='--')
ax2.set_xlabel('Angle (degrees)')
ax2.set_ylabel('Residual')
ax2.legend()
plt.tight_layout()
plt.savefig('compton_scattering.png')
