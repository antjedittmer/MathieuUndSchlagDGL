MathieuStruttscheKarte2D.m (originally BerechnungMathieu.m; combines functions BerechnungMathieu_ADshort.m and correctImagValues.m.)
Calculates real and imaginary part of characteristic exponents
Outputs excel file with characteristic mulitpliers and exponents at increasing amplification of constant (nu_0^2) and periodic (nu_c^2) part
Outputs plots STRUTTscheKarte_D*_SW*_unt*.png (Characteristic exponents: area of stability, real and imaginary part) and CharExp_D*_SW*.png (char. Exponents: Mathieu equation in title, real and imaginary part)
 
MathieuStabilitaetsgebirge3D.m (originally BerechnungMathieu.m; contains script BerechnungMathieuReSR_AD.m)
Calculates real part of characteristic exponents corresponding to different combination of nu_0^2 and nu_c^2
Outputs plots Stabilitaetsgebirg*Real*.png (Real part of char. Exponents on the z axis, nu_0^2 and nu_c^2 on the x and y axes )
 
MathieuDGL.m 
Implements the Mathieu ODE
 
MathieuTransitionMatrix.m (originally MathieuDiagramme.m; contains MathieuDiagrammeTransitionsmatrizen_GeschwX1PosX2.m as a script.
Solves the Mathieu differential equations for different damping parameters (D) and computes the evolution of the transition (monodromy) matrices over one period, generating 2D and 3D visualizations of the resulting state trajectories.
Outputs 3D (and optionally 2D) plots showing the relationship between the phase variable (), the state (), and its derivative (')

Frequency_Peters.m (Floquet Frequency Plot for Mathieu's Equation following Peters' convention)
Calculates Floquet frequencies (imaginary parts of characteristic exponents) for the Mathieu equation and separates them into branches labeled by integer multiple 
m∈[−4,4]. Outputs data files storing ϵ values and normalized frequencies ω/Ω for each branch m (Excel .xlsx and .mat formats).
​Outputs plots PetersFrequency_wXpX.png (Frequency vs. ϵ, colored lines with branch labels m, annotations for branch groups like [-1/+0], [-2/+1]).

PetersHarmonicParticipation.m (Computes normalized harmonic modal participation vs epsilon for w = 0.3/0.5/0.7 using mode tracking/continuation)
Calculates harmonic modal participation coefficients for Floquet modes of Mathieu equation x¨+(w^2+ϵsin(Ωt))x=0.
Uses mode tracking (continuation): starts from initial Floquet exponent η_0 =−wΩ, solves monodromy matrix Φ(T) via ode45 for each ϵ, selects closest eigenvalue/mode to previous step, extracts periodic mode shape Q(t), computes FFT spectrum over one period, normalizes harmonic magnitudes at frequencies 
m/T for m∈[−3,3] to get participation.
Outputs data files storing ϵ and normalized participation coefficients ϕ_m for each harmonic branch m (Excel .xlsx and .mat formats in dataFolder/Peters_HarmonicParticipation_w_XpX.*).
Outputs plots PetersHarmonicparticipatio_wXpX.png (Modal participation vs.ϵ, colored lines for each m, NaN-inserted at jumps >0.2, branch labels like m=+1 → ω/Ω ≈ 1.3, annotations marking mode switches like [-1/+0], [-2/+1]).

PetersModalWaveforms.m (Modal waveforms Q(t) for specific (w,ε) parameter combinations)
Computes periodic Floquet modal waveforms Q(t) = x(t)exp(-σt) for Mathieu equation at 6 test cases: (w,ε) = (0.7,0), (0.7,1), (0.7,2), (0.7,3), (1.0,0), (1.0,3.5).
Method: Computes monodromy matrix Φ(T) by integrating 2 column solutions from initial conditions [1;0] and [0;1], extracts Floquet exponents η and eigenvectors V. Selects stable mode, integrates over t∈[0,4π] to get x(t), extracts Q(t).
Special Peters cases (w=0.7,ε=3.0 and w=1.0,ε=3.5): Uses Peters' mode-tracking approach - solves state-transition matrix equation dΦ/dt = D(t)Φ with Φ(0)=I, extracts correct stable Floquet mode by closest match to η_target = -wΩ, constructs Q(t) = Φ(t)v_mode exp(-η_mode t) over one period, tiles to cover t∈[0,4π].
Outputs plots PetersModelWaveForm_wXpX_epsYpY.png (real(Q(t)) over t∈[0,12.5], 6 colored curves, high precision ode45 with RelTol=1e-12).

PetersRotorFlapping.m (Damped Flapping Frequency Plot - Peters' full flapping coefficients)
Method: Direct Floquet analysis of rotor blade flapping equation with Peters' coefficients C(t), K(t); computes monodromy matrix Φ(T) via ode45 on state transition equation, extracts η = log(eig(Φ(T)))/T, separates Im(η)/Ω into branches m∈[-4,4] using principal folding.
Parameters: p=1.0, γ=12.0,Ω=1 (1/rev), sweeps advance ratio μ∈[0,3].
Separates frequencies into branches m ∈ [-4,4] using principal value folding basis_freq = mod(ω/Ω + 0.5, 1) - 0.5, then ω_m/Ω = basis_freq + m.
​Outputs plots Petersrotor–bladeflapping<Style>.png (Frequency vs. μ, colored/black dots for branches [-4]...[+4], annotations [-4/+3], [-3/+2], [-1/+0], secondary [-4/+4] etc.)

compareImagCharExpFreqParticipation.m (Peters vs Arnold characteristic exponents + participation - D=0.15)
Loads Arnold reference data STRUTTscheKarte_D1dot5e-01_SW1dot0e-01_unt0.mat, computes Peters Floquet analysis for damped Mathieu ẍ + 2Dẋ + (ν+νcos(t))x = 0 with mode tracking from initial η=√ν.
Calculates composite frequency ω=sum(ω_m φ_m), growth rate σ=Re(η), harmonic participation φ_m∈[-4,4] via FFT of periodic mode shape Q(t)=Φ(t)v_mode exp(-η t).
Outputs plots:real_vs_imaginary_exponents_Mathieu_Arnold.png (σ vs ω scatter, Arnold reference)
compare_char_exp_freq_participation_w1dot0.png (4-panel: Peters/Arnold ω+σ comparison, Δω, branch frequencies m=-4:4, participation φ_m)
real_vs_imaginary_exponents_Mathieu_Peters.png (σ vs ω scatter, Peters method)

comparison_additionfactor_m.m (Arnold vs Peters harmonic frequencies)
Loads Arnold reference STRUTTscheKarte_D1dot5e-01_SW1dot0e-01_unt0.mat, auto-detects Peters harmonic data from .mat files or manual input, interpolates to common grid ν∈[0,max].
Computes: Peters-Arnold differences Δω, mean offset M≈0.XXX, stability bubbles Re(s_R), corrected Peters ω-M vs Arnold.
Outputs 4-panel plot (bubble behavior, Δω→M, Arnold stability Re(s_R), Peters corrected), reveals Peters frequency wiggles in instability bubbles while Arnold stays flat.

compareImagcharExp_PetersandArnolds.m (Harmonic analysis for undamped Mathieu w=0.3)
Method: Floquet mode tracking for ẍ + (ε+εcos(t))x = 0 (D=0), sweeps ε∈ (500 points), computes composite frequency ω=∑ω_m φ_m, growth rate Re(η), harmonic participation φ_m∈[-5,5] via FFT of periodic mode Q(t)=Φ(t)v_mode exp(-η t).
​Outputs plots (3-panel: composite freq+growth rate, area/line harmonic participation m=-5:5)
