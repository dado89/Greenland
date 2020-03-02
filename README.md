# The project in short
Implementation of Gaussian Process Emulation to reconstruct the ice morphology (global volume and displacement) of Greenland during the warmer-than-today Last Interglacial, around 120k years ago.

Collaboration with Dr Louise Sime - British Antarctic Surey, Cambridge (UK).


# Slides
Slide on the topic provided (talk given at internal PhD seminar in Leeds).

====================================================================================

Also provided a picture (Posterior_Density_3x3.pdf) comparing prior and posterior densities for the data-compatible morphologies.


# Code
A detailed description of each piece of code is given within the latter. Here a summary:
- pca_greenland.m: extracts 13 Principal Components (PCs) highlighting the regions of larger morphology variability;
- build_shapes.m:  constructs new morphologies as affine combination of the PCs;
- emul.m:          implements emulation on the dataset. That is, it predicts the simulator output (delta^18 O) for a number of morphologies, at a rate around 10^9 times faster (miliseconds rather than weeks) than the time needed to actually run the simulator;
- cross_val.m:     needed to estimated the goodness of the built emulator, with fixed correlation lengths and observational variance;
- data_match.m:    based on the emulator prediction and ice-core records, selects data-compatible Greenland Ice Sheet morphologies;
