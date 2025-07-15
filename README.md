# Compton Experiment

This repository contains the data, analysis code, and results for a Compton scattering experiment. It includes the characterization of the electronic chain, data analysis scripts, simulations, and final reports.

## Repository Structure

Here is an overview of the main directories and their contents:

-   **.vscode/**: Visual Studio Code workspace settings.
-   **`Charatterization_elettronic_chain/`**: Contains all data and analysis related to the characterization of the experimental setup's electronic components. This includes studies on the ADC, coincidence timing, efficiency, HV & gain calibration, and shaping time.
-   **`Codes/`**: Houses all the Python source code for the analysis.
    -   `setup.py`: A setup script to install the project's custom library (`libraries_x_Compton`).
    -   `ADC_studies/`: Scripts for analyzing ADC performance.
    -   `auto_coincidence/`: Scripts for automated coincidence event analysis.
    -   `Christmas_analysis/`: Specific analysis scripts, likely from a particular measurement run.
    -   `data_analysis/`: General-purpose data analysis scripts.
    -   `HV_gain_analysis/`: Scripts for high-voltage and gain calibration analysis.
-   **`Measurments/`**: Contains the raw files and plots from the experimental runs.
-   **`Pdf/`**: Stores relevant PDF documents.
-   **`Results/`**: Contains the final, processed results derived from the data analysis.
-   **`Simulations/`**: Contains simulation code, scripts, and resulting data used to model the experiment.
-   **`Target_structure/`**: Information, diagrams, or data related to the physical structure of the experimental target.
