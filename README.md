# LCMS_sim Course Package
This project simulates a bottom-up proteomics workflow, covering the main steps from protein digestion to liquid chromatography (LC) and mass spectrometry (MS1/MS2).

## It includes:
- Reading protein sequences from FASTA files
- Protein digestion into peptide fragments
- Prediction of peptide LC retention times
- Calculation of peptide molecular masses and mass-to-charge ratios (m/z)
- Generation of b- and y-ion fragments for MS/MS analysis
- Visualization of simulated lc retention times and mass spectra

## Requirements
- numpy
- pandas
- matplotlib
- pyteomics
- ipykernel
- pytest

## Installation steps
Clone the repository and install the requirements
pip install -r requirements.txt
pip install .

## HOW TO
Provide a FASTA file and use it as input to run the provided functions and simulate the complete bottom-up LC-MS workflow, from protein digestion to peptide analysis and MS1/MS2 spectrum visualization. 
The tutorials folder gives a short example.
