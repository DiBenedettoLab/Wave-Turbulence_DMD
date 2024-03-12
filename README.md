# Ocean Wave-Turbulence Decomposition MATLAB Script

## References
See [manuscript](https://arxiv.org/abs/2403.00223) for a detailed description of the method and how to use it.

```
Wave-Turbulence Decomposition MATLAB Script
│
├── Dependencies:
│   ├── MATLAB Signal Processing Toolbox 23.2 or later
│   │
│   ├── HankelMatrix.m
│   ├── DMD_processing.m
│   │   ├── DMDselective.m
│   │   └── HankelMatrix.m
│   ├── DMD_spectrum.m
│   ├── find_in_range.m
│   ├── DMD_reconstruct.m
│   │   └── textprogressbar.m
│   └── eigvals_lambda_omega.m
│
└── Main Script:
    ├── Load data (SynthData.mat)
    ├── Adjust parameters (nfft, win, fs, waveRange)
    ├── Visualize power spectra
    │
    ├── Proper Orthogonal Decomposition (POD)
    │   └── HankelMatrix.m
    ├── Dynamic Mode Decomposition (DMD)
    │   ├── DMD_processing.m
    │   │   ├── DMDselective.m
    │   │   └── HankelMatrix.m
    │   └── DMD_spectrum.m
    │
    ├── Check energy distribution
    │   └── eigvals_lambda_omega.m
    │
    └── Plotting
```

## What does this script do?
This MATLAB script helps you separate wave and turbulence motions from raw flow measurement data using dynamic mode decomposition (DMD).

## Input:
- **Raw Data:** Load your flow measurement data.
- **Sampling Frequency (fs):** Specify the frequency at which measurements were taken.
- **Wave Frequency Range:** Define the range of frequencies where waves are expected in a power density spectra of the raw data.
- **Spectral Parameters:** Set the number of FFT points (``nfft``) and window size (``win``) for spectral visualization.

## Output:
- **Out:** Structure that contains the outputs from DMD analysis such as the eigen modes (``Out.Phi``), discrete and continuous eigenvalues (``Out.lambda`` and ``Out.omega``), amplitudes (``Out.b``), and DMD separation of the raw data (``Out.Xdmd``).

## How to use:

1. **Load Your Data:** Load your flow measurement data as a $1 \times N$ array.
2. **Adjust Parameters:** Tune parameters like sampling frequency, wave frequency range, time-delay, and spectral parameters according to your data characteristics.
3. **Get the POD spectra:** The POD spectra helps determine the rank truncation of the time-delayed matrix.
4. **Run DMD:** Decompose the time-delayed matrix using DMD and reconstruct the time series in time domain.

## Custom functions
- **``HankelMatrix.m``:** Constructs Hankel matrix
- **``DMD_processing.m``:** Processes raw signal and performs DMD
- **``DMD_spectrum.m``:** Computes DMD spectrum
- **``find_in_range.m``:** Finds indexes within a specified range
- **``DMD_reconstruct.m``:** Reconstructs signals using DMD
- **``eigvals_lambda_omega.m``:** Extracts components from eigenvalues for visualization
