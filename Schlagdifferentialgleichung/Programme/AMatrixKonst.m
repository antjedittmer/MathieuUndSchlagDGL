%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Eliminiert die Periodizitäten durch symbolische Integration und erzeugt
% die A-Matrizen mit konstanten Koeffizienten
% Schreibt diese in WorkspaceA.mat
%
% benötigtes Addon in MATLAB: Symbolic MathToolbox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

syms gamma d4 d3 ebeta mu nu0 d2 psi
pi = 3.14159265358979323846;

%% 3 Blatt Rotor
    dm11 = (gamma*(d4+d3*ebeta))/2;
    dm12 = 0;
    dm13 = (gamma*mu*d3)/4;
    dm21 = 0;
    dm22 = (gamma*mu*d3*sin(3*psi)+2*gamma*(d4+d3*ebeta))/4;
    dm23 = (-gamma*mu*d3*cos(3*psi)+8)/4;
    dm31 = (gamma*mu*d3)/2;
    dm32 = (-gamma*mu*d3*cos(3*psi)-8)/4;
    dm33 = (-gamma*mu*d3*sin(3*psi)+2*gamma*(d4+d3*ebeta))/4;
    
    km11 = nu0^2;
    km12 = (gamma*mu^2*d2*sin(3*psi)+2*gamma*mu*d2*ebeta)/8;
    km13 = (-gamma*mu^2*d2*cos(3*psi))/8;
    km21 = (gamma*mu^2*d2*sin(3*psi)+2*gamma*mu*(d3+d2*ebeta))/4;
    km22 = (gamma*mu*cos(3*psi)*(2*d3+d2*ebeta)+4*(nu0^2-1))/4;
    km23 = (gamma*mu^2*d2+2*gamma*mu*sin(3*psi)*(2*d3+d2*ebeta)+4*gamma*(d4+d3*ebeta))/8;
    km31 = (-gamma*mu^2*d2*cos(3*psi))/4;
    km32 = (gamma*mu^2*d2+2*gamma*mu*sin(3*psi)*(2*d3+d2*ebeta)-4*gamma*(d4+d3*ebeta))/8;
    km33 = (gamma*mu*cos(3*psi)*(-2*d3-d2*ebeta)+4*(nu0^2-1))/4;
    
    A = -[dm11,dm12,dm13,km11,km12,km13;
          dm21,dm22,dm23,km21,km22,km23;
          dm31,dm32,dm33,km31,km32,km33;
          -1,0,0,0,0,0;
          0,-1,0,0,0,0;
          0,0,-1,0,0,0];

mat = A;

sympref('FloatingPointOutput',true);

A3 = 1/(2*pi) * int(mat,psi,0,2*pi)
subs A3


%% 4 Blatt Rotor
    dm11 = (gamma*(d3*ebeta+d4))/2;
    dm12 = 0;
    dm13 = (d3*gamma*mu)/4;
    dm14 = 0;
    dm21 = 0;
    dm22 = (gamma*(d3*ebeta+d4))/2;
    dm23 = 2;
    dm24 = (-sin(2*psi)*d3*gamma*mu)/2;
    dm31 = (d3*gamma*mu)/2;
    dm32 = -2;
    dm33 = (gamma*(d3*ebeta+d4))/2;
    dm34 = (cos(2*psi)*d3*gamma*mu)/2;
    dm41 = 0;
    dm42 = (-sin(2*psi)*d3*gamma*mu)/4;
    dm43 = (cos(2*psi)*d3*gamma*mu)/4;
    dm44 = (gamma*(d3*ebeta+d4))/2;
    
    km11 = nu0^2;
    km12 = (d2*ebeta*gamma*mu)/4;
    km13 = 0;
    km14 = (-sin(2*psi)*d2*gamma*mu^2)/4;
    km21 = (gamma*mu*(d2*ebeta+d3))/2;
    km22 = (sin(4*psi)*d2*gamma*mu^2+8*nu0^2-8)/8;
    km23 = (gamma*(-cos(4*psi)*d2*mu^2+d2*mu^2+4*d3*ebeta+4*d4))/8;
    km24 = (-cos(2*psi)*gamma*mu*(d2*ebeta+d3))/2;
    km31 = 0;
    km32 = (gamma*(-cos(4*psi)*d2*mu^2+d2*mu^2-4*d3*ebeta-4*d4))/8;
    km33 = (-sin(4*psi)*d2*gamma*mu^2+8*nu0^2-8)/8;
    km34 = (-sin(2*psi)*gamma*mu*(d2*ebeta+d3))/2;
    km41 = (-sin(2*psi)*d2*gamma*mu^2)/4;
    km42 = (cos(2*psi)*gamma*mu*(-d2*ebeta-2*d3))/4;
    km43 = (sin(2*psi)*gamma*mu*(-d2*ebeta-2*d3))/4;
    km44 = nu0^2;
    
    A = -[dm11,dm12,dm13,dm14,km11,km12,km13,km14;
          dm21,dm22,dm23,dm24,km21,km22,km23,km24;
          dm31,dm32,dm33,dm34,km31,km32,km33,km34;
          dm41,dm42,dm43,dm44,km41,km42,km43,km44;
          -1,0,0,0,0,0,0,0;
          0,-1,0,0,0,0,0,0;
          0,0,-1,0,0,0,0,0;
          0,0,0,-1,0,0,0,0];

mat = A;

sympref('FloatingPointOutput',true);

A4 = 1/(2*pi) * int(mat,psi,0,2*pi)
subs A4

%% 5 Blatt Rotor
    dm11 = (gamma*(d4+d3*ebeta))/2;
    dm12 = 0;
    dm13 = (gamma*mu*d3)/4;
    dm14 = 0;
    dm15 = 0;
    dm21 = 0;
    dm22 = (gamma*(d4+d3*ebeta))/2;
    dm23 = 2;
    dm24 = 0;
    dm25 = (gamma*mu*d3)/4;
    dm31 = (gamma*mu*d3)/2;
    dm32 = -2;
    dm33 = (gamma*(d4+d3*ebeta))/2;
    dm34 = (-gamma*mu*d3)/4;
    dm35 = 0;
    dm41 = 0;
    dm42 = 0;
    dm43 = (-gamma*mu*d3)/4;
    dm44 = (gamma*mu*d3*sin(5*psi)+2*gamma*(d4+d3*ebeta))/4;
    dm45 = (-gamma*mu*d3*cos(5*psi)+16)/4;
    dm51 = 0;
    dm52 = (gamma*mu*d3)/4;
    dm53 = 0;
    dm54 = (-gamma*mu*d3*cos(5*psi)-16)/4;
    dm55 = (-gamma*mu*d3*sin(5*psi)+2*gamma*(d4+d3*ebeta))/4;
    
    km11 = nu0^2;
    km12 = (gamma*mu*d2*ebeta)/4;
    km13 = 0;
    km14 = 0;
    km15 = (gamma*mu^2*d2)/8;
    km21 = (gamma*mu*(d3+d2*ebeta))/2;
    km22 = nu0^2-1;
    km23 = (gamma*mu^2*d2+4*gamma*(d4+d3*ebeta))/8;
    km24 = (gamma*mu^2*d2*sin(5*psi)+2*gamma*mu*(-d3+d2*ebeta))/8;
    km25 = (-gamma*mu^2*d2*cos(5*psi))/8;
    km31 = 0;
    km32 = (gamma*mu^2*d2-4*gamma*(d4+d3*ebeta))/8;
    km33 = nu0^2-1;
    km34 = (-gamma*mu^2*d2*cos(5*psi))/8;
    km35 = (-gamma*mu^2*d2*sin(5*psi)+2*gamma*mu*(-d3+d2*ebeta))/8;
    km41 = 0;
    km42 = (gamma*mu^2*d2*sin(5*psi)+2*gamma*mu*(2*d3+d2*ebeta))/8;
    km43 = (-gamma*mu^2*d2*cos(5*psi))/8;
    km44 = (gamma*mu*cos(5*psi)*(3*d3+d2*ebeta)+4*(nu0^2-4))/4;
    km45 = (gamma*mu*sin(5*psi)*(3*d3+d2*ebeta)+4*gamma*(d4+d3*ebeta))/4;
    km51 = (gamma*mu^2*d2)/4;
    km52 = (-gamma*mu^2*d2*cos(5*psi))/8;
    km53 = (-gamma*mu^2*d2*sin(5*psi)+2*gamma*mu*(2*d3+d2*ebeta))/8;
    km54 = (gamma*mu*sin(5*psi)*(3*d3+d2*ebeta)-4*gamma*(d4+d3*ebeta))/4;
    km55 = (gamma*mu*cos(5*psi)*(-3*d3-d2*ebeta)+4*(nu0^2-4))/4;
    
    A = -[dm11,dm12,dm13,dm14,dm15,km11,km12,km13,km14,km15;
          dm21,dm22,dm23,dm24,dm25,km21,km22,km23,km24,km25;
          dm31,dm32,dm33,dm34,dm35,km31,km32,km33,km34,km35;
          dm41,dm42,dm43,dm44,dm45,km41,km42,km43,km44,km45;
          dm51,dm52,dm53,dm54,dm55,km51,km52,km53,km54,km55;
          -1,0,0,0,0,0,0,0,0,0;
          0,-1,0,0,0,0,0,0,0,0;
          0,0,-1,0,0,0,0,0,0,0;
          0,0,0,-1,0,0,0,0,0,0;
          0,0,0,0,-1,0,0,0,0,0];

mat = A;

sympref('FloatingPointOutput',true);

A5 = 1/(2*pi) * int(mat,psi,0,2*pi)
subs A5

save("WorkspaceA.mat","A3","A4","A5")


