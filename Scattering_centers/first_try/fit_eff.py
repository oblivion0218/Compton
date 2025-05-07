import numpy as np
import matplotlib.pyplot as plt
from iminuit import Minuit
from iminuit.cost import LeastSquares

# Data provided by the user
angles_deg = np.array([15, 35, 40, 50, 60, 70, 90, 110])
efficiencies = np.array([0.34967, 0.40536, 0.42514, 0.47027, 0.52115, 0.57572, 0.68808, 0.79343])

# Convert angles to radians
angles_rad = np.radians(angles_deg)

# Define the exponential function to fit
# f(x) = a * exp(-b * x) + c
def model_func(x, a, b, c):
    return a * np.exp(-b * x) + c

# Create the least squares cost function
# We need error estimates for a proper fit, but let's assume equal errors for now
# If you have errors (e.g., y_err), pass them to LeastSquares:
# least_squares = LeastSquares(angles_rad, efficiencies, y_err, model_func)
least_squares = LeastSquares(angles_rad, efficiencies, 1.0, model_func) # Assuming error=1 for all points

# Instantiate Minuit
# Provide initial guesses for parameters a, b, c
# Guesses: a ~ range/2, b ~ 1/range_x, c ~ min_y
m = Minuit(least_squares, a=0.5, b=1.0, c=0.3)

# Run the minimization
m.migrad()

# Check the fit status
m.hesse() # Calculate uncertainties accurately
print(m.fmin) # Print fit status information
print(m.params) # Print fitted parameters and uncertainties

# Plot the results
plt.figure(figsize=(8, 6))
plt.scatter(angles_rad, efficiencies, label='Data points', color='blue')

# Generate points for the fitted curve
x_fit = np.linspace(min(angles_rad), max(angles_rad), 200)
y_fit = model_func(x_fit, *m.values) # Use fitted parameter values

plt.plot(x_fit, y_fit, label=f'Fit: a*exp(-b*x)+c\na={m.values["a"]:.3f}±{m.errors["a"]:.3f}\nb={m.values["b"]:.3f}±{m.errors["b"]:.3f}\nc={m.values["c"]:.3f}±{m.errors["c"]:.3f}', color='red')

plt.xlabel("Angle (radians)")
plt.ylabel("Efficiency")
plt.title("Fit of Efficiency vs Angle")
plt.legend()
plt.grid(True)
plt.savefig("efficiency_fit.png")

# Print fitted parameters again for clarity
print("\nFitted Parameters:")
for p in m.parameters:
    print(f"{p} = {m.values[p]:.5f} +/- {m.errors[p]:.5f}")