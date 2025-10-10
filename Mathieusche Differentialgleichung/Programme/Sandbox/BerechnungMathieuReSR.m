%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Erzeugen des Stabilitaetsgebirges in den Grenzen von nu_02 und nu_C2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc ; clear all ; close all;

%Schrittweite
SW = 0.1;

%Grenzen fuer nu_02 und nu_C2
unt0 = 0;
untC = 0;
ob0 = 20;
obC = 20;

%Zaehlveriablen
o=1;
AnzGl=2;

%Parameter
D = 0.001;
t0 = 0.0;
T = 2*pi;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Vorbereitung der in den Schleifen zu fuellenden Arrays zwecks
%Programmbeschleunigung
Diagonal=diag(ones(AnzGl,1));
Monodromie = zeros(AnzGl);
CharEx = zeros(length(unt0:SW:ob0),AnzGl*3);
plotwert = zeros(length(unt0:SW:ob0)*length(untC:SW:obC),3);
Z=zeros(length(unt0:SW:ob0),length(untC:SW:obC));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[X,Y] = meshgrid(unt0:SW:ob0,untC:SW:obC);

for m = unt0:SW:ob0
    l=1;
    for n = untC:SW:obC   

        nu_02 = m;
        nu_C2 = n;
        
        %Optionen fuer die Genauigkeit und Toleranz fuer den ODE-Solver
        options = odeset('RelTol',1e-10,'AbsTol',1e-12);

        for k=1:AnzGl
            sol = ode45(@(psi,x)MathieuDGL(psi,x,D,nu_02,nu_C2),[t0,T],Diagonal(:,k),options);
            MonoVek = deval(sol,T);
            Monodromie(:,k) = MonoVek;
        end
        
        %charakteristische Multiplikatoren (Eigenwerte der Monodromiematrix)
        p = [1 -(Monodromie(1,1)+Monodromie(2,2)) (Monodromie(1,1)*Monodromie(2,2))-Monodromie(1,2)*Monodromie(2,1)];
        e = roots(p).';
        
        ReLambda = 1/T * log(abs(e(1,1)));

        Z(l,o)= ReLambda;
        l=l+1;  
    end
   o=o+1;
end

%% Grafische Darstellung

surf(X,Y,Z)
view([5 14])
xlim([0 20]);
ylim([0 20]);
zlim([-2 2]);
xlabel('$\rm{Parameter} \; \nu_0^2 \rm{[-]}$','interpreter','latex','FontSize', 24);
ylabel('$\rm{Parameter} \; \nu_C^2 \rm{[-]}$','interpreter','latex','FontSize', 24);
zlabel('$Re(s_R) \rm{[-]}$','interpreter','latex','FontSize', 24);

