I have compared the following 2 approaches for the computation of the 95% cl L.L.:

 - Approach1: Consider only the cases where the cutoff is not detected. In this case, I re-normalized the distribution of non-detected cutoffs (and all other parameters) to 1000, by simply running more simulations and selecting the cases of non-detection.
 - Approach2: Consider the full distribution of 1000 fitted cutoffs, thus including both detected and undetected cutoffs.
 
In this folder, you see the results for a bright source (50 mCrab, -1.7, 200 TeV), with size 0.1 deg.
