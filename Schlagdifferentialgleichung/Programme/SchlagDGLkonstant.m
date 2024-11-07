%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Funktion mit dem zugrundeliegenden konstanten Schlag-
% Differentialgleichungssystem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xpunkt,A] = SchlagDGLkonstant(~,x,gamma,d2,d3,d4,mu,ebeta,nu0,Blatt)

if Blatt == 3
    dm11 = -0.5000*gamma*(d4 + d3*ebeta);
    dm12 = 0;
    dm13 = -0.2500*d3*gamma*mu;
    dm21 = 0;
    dm22 = -0.5000*gamma*(d4 + d3*ebeta);
    dm23 = -2.0000;
    dm31 = -0.5000*d3*gamma*mu;
    dm32 = 2.0000;
    dm33 = -0.5000*gamma*(d4 + d3*ebeta);
    
    km11 = -1.0000*nu0^2;
    km12 = -0.2500*d2*ebeta*gamma*mu;
    km13 = 0;
    km21 = -0.5000*gamma*mu*(d3 + d2*ebeta);
    km22 = 1.0000 - 1.0000*nu0^2;
    km23 = -0.1250*gamma*(d2*mu^2 + 4*d4 + 4*d3*ebeta);
    km31 = 0;
    km32 = 0.1250*gamma*(- d2*mu^2 + 4*d4 + 4*d3*ebeta);
    km33 = 1.0000 - 1.0000*nu0^2;
    
    A = -[dm11,dm12,dm13,km11,km12,km13;
          dm21,dm22,dm23,km21,km22,km23;
          dm31,dm32,dm33,km31,km32,km33;
          -1,0,0,0,0,0;
          0,-1,0,0,0,0;
          0,0,-1,0,0,0];
elseif Blatt == 4
    dm11 = -0.5000*gamma*(d4 + d3*ebeta);
    dm12 = 0;
    dm13 = -0.2500*d3*gamma*mu;
    dm14 = 0;
    dm21 = 0;
    dm22 = -0.5000*gamma*(d4 + d3*ebeta);
    dm23 = -2;
    dm24 = 0;
    dm31 = -0.5000*d3*gamma*mu;
    dm32 = 2;
    dm33 = -0.5000*gamma*(d4 + d3*ebeta);
    dm34 = 0;
    dm41 = 0;
    dm42 = 0;
    dm43 = 0;
    dm44 = -0.5000*gamma*(d4 + d3*ebeta);
    
    km11 = -1.0000*nu0^2;
    km12 = -0.2500*d2*ebeta*gamma*mu;
    km13 = 0;
    km14 = 0;
    km21 = -0.5000*gamma*mu*(d3 + d2*ebeta);
    km22 = 1.0000 - 1.0000*nu0^2;
    km23 = -0.1250*gamma*(d2*mu^2 + 4*d4 + 4*d3*ebeta);
    km24 = 0;
    km31 = 0;
    km32 = 0.1250*gamma*(- d2*mu^2 + 4*d4 + 4*d3*ebeta);
    km33 = 1.0000 - 1.0000*nu0^2;
    km34 = 0;
    km41 = 0;
    km42 = 0;
    km43 = 0;
    km44 = -1.0000*nu0^2;
    
    A = [dm11,dm12,dm13,dm14,km11,km12,km13,km14;
          dm21,dm22,dm23,dm24,km21,km22,km23,km24;
          dm31,dm32,dm33,dm34,km31,km32,km33,km34;
          dm41,dm42,dm43,dm44,km41,km42,km43,km44;
          1,0,0,0,0,0,0,0;
          0,1,0,0,0,0,0,0;
          0,0,1,0,0,0,0,0;
          0,0,0,1,0,0,0,0];
elseif Blatt == 5
    dm11 = -0.5000*gamma*(d4 + d3*ebeta);
    dm12 = 0;
    dm13 = -0.2500*d3*gamma*mu;
    dm14 = 0;
    dm15 = 0;
    dm21 = 0;
    dm22 = -0.5000*gamma*(d4 + d3*ebeta);
    dm23 = -2.0000;
    dm24 = 0;
    dm25 = -0.2500*d3*gamma*mu;
    dm31 = -0.5000*d3*gamma*mu;
    dm32 = 2.0000;
    dm33 = -0.5000*gamma*(d4 + d3*ebeta);
    dm34 = 0.2500*d3*gamma*mu;
    dm35 = 0;
    dm41 = 0;
    dm42 = 0;
    dm43 =  0.2500*d3*gamma*mu;
    dm44 = -0.5000*gamma*(d4 + d3*ebeta);
    dm45 = -4.0000;
    dm51 = 0;
    dm52 = -0.2500*d3*gamma*mu;
    dm53 = 0;
    dm54 = 4.0000;
    dm55 = -0.5000*gamma*(d4 + d3*ebeta);
    
    km11 = -1.0000*nu0^2;
    km12 = -0.2500*d2*ebeta*gamma*mu;
    km13 = 0;
    km14 = 0;
    km15 = -0.1250*d2*gamma*mu^2;
    km21 = -0.5000*gamma*mu*(d3 + d2*ebeta);
    km22 = 1.0000 - 1.0000*nu0^2;
    km23 = -0.1250*gamma*(d2*mu^2 + 4*d4 + 4*d3*ebeta);
    km24 = 0.2500*gamma*mu*(d3 - d2*ebeta);
    km25 = 0;
    km31 = 0;
    km32 = 0.1250*gamma*(- d2*mu^2 + 4*d4 + 4*d3*ebeta);
    km33 = 1.0000 - 1.0000*nu0^2;
    km34 = 0;
    km35 = 0.2500*gamma*mu*(d3 - d2*ebeta);
    km41 = 0;
    km42 = -0.2500*gamma*mu*(2*d3 + d2*ebeta);
    km43 = 0;
    km44 = 4.0000 - 1.0000*nu0^2;
    km45 = -1.0000*gamma*(d4 + d3*ebeta);
    km51 = -0.2500*d2*gamma*mu^2;
    km52 = 0;
    km53 = -0.2500*gamma*mu*(2*d3 + d2*ebeta);
    km54 = 1.0000*gamma*(d4 + d3*ebeta);
    km55 = 4.0000 - 1.0000*nu0^2;
    
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
end


xpunkt = A*x;
end
