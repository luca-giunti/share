Each folder corresponds to a different set of injected source parameters. The naming scheme is:

        {norm at 1 TeV [mCrab]}_{spectral index}_{cutoff energy [TeV]}.
        
Each folder contains 1 text file (reporting the distribution of TS) and 2 subfolders: 
    - "PL" contains the distribution of fitted parameters for the Power Law model (i.e. `norm`, `index`, `RA`, `DEC`, `bkg_norm`, `bkg_tilt`) 
    - "ECPL" contains the distribution of fitted parameters for the Ecponential Cutoff Power Law model (i.e. `norm`, `index`, `lambda_`, `RA`, `DEC`, `bkg_norm`, `bkg_tilt`) 

The injected `RA` and `DEC` values are, respectively: 75.048 deg, 60.075 deg. However, this is not very relevant, and in the text files I just wrote the difference between the best-fit `RA`/`DEC` values and the injected ones: `RA_fitted - RA_simulated`, `DEC_fitted - DEC_simulated` .

**Regarding the TS computation**: fit statistics in gammapy are all log-likelihood functions normalized like chi-squares, i.e. if `L` is the likelihood function used, they follow the expression `- 2 log L`. In our specific case, the fit likelihood is computed using Cash statistic (reference: http://cxc.harvard.edu/sherpa/ahelp/cash.html), and corresponds to the following TS definition:

        TS = -2 (log(like_PL) - log(like_ECPL)) = Cash_PL - Cash_ECPL 