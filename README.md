# Study of Compton Effect with NaI Detector

This repository contains the code and images related to our work in a particle physics laboratory, where we are studying the **Compton Effect** using two **NaI scintillation detectors** and a radioactive source of **Na-22**.

## Overview

The purpose of this project is to explore and analyze the Compton scattering of gamma rays, focusing on the following aspects:
- Testing and calibrating the electronic chain.
- Determining the **best energy resolution** of the detectors (named **Franco** and **Ugo**, with Franco being the one with better resolution).
- Verifying the **linearity of the ADC** (Analog-to-Digital Converter).
- Understanding the interaction of 511 keV and 1274 keV gamma rays in the detector.

## Experimental Setup
<div align="center">
  <img src="experimental_setup.png" alt="Experimental setup" width="400">
</div>


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
