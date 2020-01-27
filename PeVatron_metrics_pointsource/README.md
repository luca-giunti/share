This folder contains the PeVatron metrics, i.e. the 95% C.L. lower limits of the fitted cutoffs distribution for a point source. The lower limits have been computed considering 1000 simulations and selecting both detected and undetected cutoffs. We used gammapy, and assumed 10h of observation time.

The `1D` and `3D` folders both contain:
- a .png file showing the metrics
- a .txt file with the data.

Each .txt file contains 4 column.
- First column: Spectral Index of the simulated source 
- Second column: Energy cutoff of the simulated source (TeV) (negative values, e.g. -1.7) of the simulated source
- Third colun: List of Flux normalizations at 1 TeV (mCrab) of the simulated source
- Fourth column: List of 95% C.L. lower limits (TeV)
