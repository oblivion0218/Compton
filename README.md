# Study of Compton Effect with NaI Detector

This repository contains the code and images related to our work in a particle physics laboratory, where we are studying the **Compton Effect** using two **NaI scintillation detectors** and a radioactive source of **Na-22**.

---

## Experimental Setup
<div align="center">
  <img src="experimental_setup.png" alt="Experimental setup" width="600">
</div>

We named the detector used for gating **Ugo** and the one dedicated to spectroscopy **Franco**.

### Components
1. **Radioactive Source (Na-22):**
   - Emits two 511 keV gamma rays in opposite directions due to positron annihilation.

2. **NaI Detectors ("Ugo" and "Franco"):**
   - Detect gamma rays at specific scattering angles.
   - Each detector is connected to a pre-amplifier for signal amplification.

3. **Pre-amplifiers:**
   - Amplify the electrical signals produced by the NaI detectors.

4. **Amplifier:**
   - Further process the signals received from the pre-amplifiers.

5. **Amplifier TSCA:**
   - Further process the signals received from the pre-amplifiers.
   - The TSCA acts as a single-channel analyzer (SCA), selecting signals based on their amplitude, which corresponds to the energy of the detected gamma ray.
   - The TSCA generates timing signals (TTL) when a valid event is detected.

6. **NIM-TTL Converter:**
   - Converts NIM signals to TTL signals for compatibility with other electronic components.

7. **Dual Timer:**
   - Generates precise timing signals to gate the data acquisition.

8. **MultiChannel Analyzer (MCA):**
   - Records the energy spectrum of detected gamma rays and visualizes it on a computer.

9. **Target (Optional):**
   - Serves as the scattering medium for gamma rays, used to measure Compton scattering.

### Signal Flow
1. Gamma rays from the Na-22 source are detected by the NaI detectors ("Ugo" and "Franco").
2. The electrical signals from the detectors are sent to their respective pre-amplifiers through coaxial cables.
3. The amplified signals are routed to the amplifier and amplifier TSCA, where further processing occurs.
4. Signals from the amplifier TSCA (TTL) is sent to the NIM-TTL converter to ensure compatibility with the Dual Timer.
5. The Dual Timer stretches the signal to make it correctly readable by the MCA.
6. Signals from the Dual Timer (NIM) is sent to the NIM-TTL converter to ensure compatibility with the MCA.
7. The processed signals are analyzed by the MCA and displayed on a computer, showing the energy spectrum of the detected events.

---

## Characterization of the electronic chain
### Choosing the best parameters for the amplifiers
We need to determine the optimal High Voltage (HV) and Gain settings for the Ugo detector, in combination with the TSCA amplifier, as well as for the Franco detector coupled with its amplifier. This ensures the detectors operate with maximum energy resolution and signal quality.
<p align="center">
  <img src="HV&Gain_calibration/Energetic_resolution/Heatmap/heatmap_Calibrazione_Franco_AMP.png" alt="Best parameters for Ugo - AMP TSCA" width="45%">
  <img src="HV&Gain_calibration/Energetic_resolution/Heatmap/heatmap_Calibrazione_Ugo_AMP-TSCA.png" alt="Best parameters for Franco - AMP" width="45%">
</p>












