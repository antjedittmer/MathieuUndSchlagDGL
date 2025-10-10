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