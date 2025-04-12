import numpy as np
import matplotlib.pyplot as plt

# data = {angle: [[energy, energy_error], [rate, rate_error]]}
data = {40: [[411.76, 0.35], [0.0946, 0.0024]], 60: [[339.31, 0.30], [0.0452, 0.0013]], 70: [[306.51, 0.15], [0.0410, 0.0006]],
        90: [[254.53, 0.11], [0.0311, 0.0005]], 110: [[217.49, 0.16], [0.0181, 0.0005]]}

S = 175000 * 903/1000 # Bq (Only for 511 KeV)
epsilon_gate = 0.4796 # Gate efficiency
solid_angle = 0.0197 # rad
scale_factor = S * epsilon_gate * solid_angle / (4 * np.pi)
print(f'Scale factor: {scale_factor}')
# Correct scale factor about 7.5e8

angle_error = np.sqrt(1 + 1 + np.arctan(0.5/25) ** 2)

def compton_energy(angle):
    return 511 / (2 - np.cos(angle * np.pi / 180))

def klein_nishina(angle, scale_factor=scale_factor):
    alpha = 1/137
    m_e = 511
    r = compton_energy(angle)/511
    return  scale_factor * (alpha ** 2 / (2 * m_e ** 2)) * r ** 2 * (r + 1 / r - np.sin(angle * np.pi / 180) ** 2)

thetas = np.linspace(0, 180, 100)

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 10))
# Energy vs angle
for angle, values in data.items():
    energy, energy_error = values[0]
    ax1.errorbar(angle, energy, yerr=energy_error, xerr=angle_error, fmt='.', label=f'{angle}°')
ax1.plot(thetas, compton_energy(thetas), label='Compton Energy', color='blue')
ax1.set_xlabel('Angle (degrees)')
ax1.set_ylabel('Energy (keV)')
ax1.set_title('Energy vs Angle')
ax1.legend()

# Rate vs angle
for angle, values in data.items():
    rate, rate_error = values[1]
    ax2.errorbar(angle, rate, yerr=rate_error, xerr=angle_error, fmt='.', label=f'{angle}°')
ax2.plot(thetas, klein_nishina(thetas, 7.5e8), label='Klein-Nishina', color='blue')
ax2.set_xlabel('Angle (degrees)')
ax2.set_ylabel('Rate (counts/s)')
ax2.set_title('Rate vs Angle')
ax2.legend()

plt.tight_layout()
plt.savefig('energy_rate_vs_angle.png')

