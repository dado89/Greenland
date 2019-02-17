# Greenland
Implementation of Gaussian Process Emulation to reconstruct the ice morphology (global volume and displacement) of Greenland during the warmer-than-today Last Interglacial, around 120k years ago.

A detailed description of each of the pieces of code is given within the code itself. In summary:
- pca_greenland.m: extracts 13 Principal Components (PCs) highlighting the regions of larger morphology variability;
- build_shapes.m:  constructs new morphologies as affine combination of the PCs;
- emul.m:          implements emulation on the dataset. That is, it predicts the simulator output (delta^18 O) for a number of morphologies, at a rate around 10^9 times faster (miliseconds rather than weeks) than the time needed to actually run the simulator;
- cross_val.m:     needed to estimated the goodness of the build emulator, with fixed correlation lengths and observational variance;
- data_match.m:    based on the emulator prediction and ice-core records, selects data-compatible Greenland Ice Sheet morphologies;
- Posterior_Density_3x3.pdf: a picture comparing prior and posterior densities for the data-compatible morphologies.
