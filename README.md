# Study of Compton Effect with NaI Detector

This repository contains the code and images related to our work in a particle physics laboratory, where we are studying the **Compton Effect** using two **NaI scintillation detectors** and a radioactive source of **Na-22**.

## Experimental Setup
<div align="center">
  <img src="experimental_setup.png" alt="Experimental setup" width="600">
</div>

### Components
1. **Radioactive Source (Na-22):**
   - Emits two 511 keV gamma rays in opposite directions due to positron annihilation.

2. **NaI Detectors ("Ugo" and "Franco"):**
   - Detect gamma rays at specific scattering angles.
   - Each detector is connected to a pre-amplifier for signal amplification.

3. **Pre-amplifiers:**
   - Amplify the electrical signals produced by the NaI detectors.

4. **Amplifiers TSCA:**
   - Further process the signals received from the pre-amplifiers.

5. **NIM-TTL Converter:**
   - Converts NIM signals to TTL signals for compatibility with other electronic components.

6. **Dual Timer:**
   - Generates precise timing signals to gate the data acquisition.

7. **MultiChannel Analyzer (MCA):**
   - Records the energy spectrum of detected gamma rays and visualizes it on a computer.

8. **Target (Optional):**
   - Serves as the scattering medium for gamma rays, used to measure Compton scattering.

### Signal Flow
1. Gamma rays from the Na-22 source are detected by the NaI detectors ("Ugo" and "Franco").
2. The electrical signals from the detectors are sent to their respective pre-amplifiers through coaxial cables.
3. The amplified signals are routed to the amplifiers TSCA, where further processing occurs.
4. Signals from the amplifiers are sent to the NIM-TTL converter to ensure compatibility with the Dual Timer.
5. The Dual Timer provides timing signals to the MCA for data acquisition.
6. The processed signals are analyzed by the MCA and displayed on a computer, showing the energy spectrum of the detected events.


1. **Source**: Sodium-22 (Na-22), emitting gamma rays of 511 keV and 1274 keV.
2. **Detectors**: 
   - **Franco**: Detector with better energy resolution.
   - **Ugo**: Secondary detector for comparison.
3. **Electronics**:
   - two PreAMP
   - AMP
   - AMP TSCA
   - ADC
   - MCA
