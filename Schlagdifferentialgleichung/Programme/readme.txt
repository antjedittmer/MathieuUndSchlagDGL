testSkriptSpiegeln.m
Processes characteristic real parts of eigenvalues (CharExRe1, CharExRe2) loaded from Workspace.mat.
Computes smooth positive and mirrored negative eigenvalue branches to correct discontinuities between them.
Generates corrected eigenvalues (CharExRe1Cor, CharExRe2Cor) ensuring continuous transitions.
Outputs two plots: original vs. corrected eigenvalues and comparison of both corrected branches.

BerechnungSchlagDGL.m 
Computes characteristic exponents and multipliers of rotor systems using Floquet theory for different rotor configurations.
Solves system matrices derived from rotor dynamics (SchlagDGL) for increasing values of the parameter μ (from MuMin to MuMax), calculates monodromy matrices, and extracts real and imaginary parts of characteristic exponents.
Applies smoothing and correction to ensure continuous eigenvalue behavior and saves results to .mat and .xlsx files.
Generates several diagnostic plots showing eigenvalue evolution, corrected characteristic exponents, and absolute values of characteristic multipliers.

BerechnungSchlagDGL1.m
Purpose: Performs a detailed stability analysis for a 3-blade see-saw rotor by calculating characteristic exponents ($s$) and multipliers ($\Lambda$) across a swept parameter ($\mu_{\text{param}} \in [0, 10]$).
Core Motivation: To produce smooth, continuous stability diagrams by overcoming numerical issues and branch crossings inherent to Floquet analysis.
Key Functionality:
Floquet Theory: Solves the system's periodic ODEs to compute the Monodromy Matrix.
Data Correction: Implements custom logic for continuous eigenvalue tracking, phase unwrapping of the imaginary part (frequency), and a mirroring technique to correct non-physical sign flips in the real part (damping).
Output: Generates stability diagrams, plots of multipliers, and exports corrected data to .mat and .xlsx files.

