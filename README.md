# Compton Effect with NaI Detector

This repository contains the code and images related to my work in a particle physics laboratory, where I am studying the **Compton Effect** using an **NaI scintillation detector** and a radioactive source of **Na-22**.

## Overview

The purpose of this project is to explore and analyze the Compton scattering of gamma rays, focusing on the following aspects:
- Testing and calibrating the electronic chain.
- Determining the **best energy resolution** of the detectors (named **Franco** and **Ugo**, with Franco being the one with better resolution).
- Verifying the **linearity of the ADC** (Analog-to-Digital Converter).
- Understanding the interaction of 511 keV and 1274 keV gamma rays in the detector.

This repository serves as a platform to share:
- The Python/C++ programs developed for data acquisition and analysis.
- Plots and images obtained from the experiments.
- Supporting documentation for experimental setups.


## Experimental Setup

1. **Source**: Sodium-22 (Na-22), emitting gamma rays of 511 keV and 1274 keV.
2. **Detectors**: 
   - **Franco**: Detector with better energy resolution.
   - **Ugo**: Secondary detector for comparison.
3. **Electronics**:
   - Preamp, shaping amplifiers, ADC modules.
   - Focused on fine-tuning for optimal performance.

## Analysis Goals

- **Energy Resolution**: Evaluating the detectors' performance in resolving close energy peaks.
- **ADC Linearity**: Ensuring that the system's ADC maintains linearity over the energy range.
- **Gamma Interaction Study**: Observing and analyzing the effects of Compton scattering on gamma-ray spectra.

## How to Use

1. Clone the repository:
   ```bash
   git clone https://github.com/oblivion0218/Compton.git
   cd Compton

