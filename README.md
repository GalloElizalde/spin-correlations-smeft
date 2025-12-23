# EFT differential fit – analysis notes & code

Repository containing an EFT / SMEFT differential analysis developed during my master thesis work. 
The purpose of this codebase is reconstruct the Spin-density-matrix of the ttbar system ant to extract 
constraints on Wilson coefficients from quantum observables.


---

## What this repo does

C++/ROOT → observables → EFT likelihood → Wilson coefficient fits → plots.

More explicitly:
- C++/ROOT macros generate histograms and quantum observables
- Julia code defines the EFT model and performs likelihood fits
- Python scripts post-process results and produce figures


---

## Structure

data/ input and intermediate datasets
ROOT/ ROOT files with the simulated data
macros/ histogram and observable construction
julia/ EFT likelihood and fitting code
python/ plotting and post-processing
plots/ final figures
results/ numerical outputs, ROOT files and intermediate figures


---
