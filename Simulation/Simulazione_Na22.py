import gamma_simulation as gs
import matplotlib.pyplot as plt

print(gs.cross_section_compton(511), gs.cross_section_photoelectron(511))


result = []
entries = 10000  # Increase number of events for better statistics

# Run detection for each gamma ray event
for _ in range(entries):
    x = gs.gamma_detection()
    if x is not None:
        for xi in x:
            result.append(xi)

# Plot histogram with a high number of bins for better spectral resolution
plt.hist(result, bins=2000, color='blue', edgecolor='blue') # range=(0, 2000)
plt.xlabel('Electron Energy (keV)')
plt.ylabel('Counts')

plt.axvline(511, color='red', linestyle='--', linewidth=1.5, label='511 keV')
plt.axvline(1274, color='red', linestyle='--', linewidth=1.5, label='1274 keV')

plt.legend()
plt.savefig("simulazione_realistica.png")  

print("\nFine\n")