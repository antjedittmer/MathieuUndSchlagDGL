%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Funktion mit dem zugrundeliegenden Differentialgleichungssystem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function A = MathieuDGLlin(D,nu_02,nu_C2)

a11 = -2*D;
a12 = -(nu_02 + nu_C2);
a21 = 1.0;
a22 = 0.0;
A = [a11, a12;a21,a22];

end