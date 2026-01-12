%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Berechnung der Struttschen Karten in den Grenzen von nu_02 und nu_C2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc ; clear all ; close all;

%Schrittweite
SW = 0.1;

%Grenzen fuer nu_02 und nu_C2
unt0 = 0;
untC = 0;
ob0 = 9;
obC = 9;

%Zaehlveriablen
l=1;
o=1;

%Anzahl der Gleichungen des DGL-Systems
Nz=2;

%Parameter
D = 0.001;
t0 = 0.0;
T = 2*pi;
tspan = t0:0.0001:T;

%Anfangsbedingungen
Diagonal=diag(ones(Nz,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Vorbereitung der in den Schleifen zu fuellenden Arrays zwecks
%Programmbeschleunigung
Monodromie = zeros(Nz);
CharMult =  zeros(length(unt0:SW:ob0),2);
CharEx = zeros(length(unt0:SW:ob0),Nz*2);
plotwertstabil = zeros(length(unt0:SW:ob0)*length(untC:SW:obC),3);
plotwertinstabil = zeros(length(unt0:SW:ob0)*length(untC:SW:obC),3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for m = unt0:SW:ob0
    for n = untC:SW:obC   

        nu_02 = m ;
        nu_C2 = n ;
        
        %Schleife intergiert Nz Mal numerisch nacheinaner fuer die
        %Spaltenvektren der Einheitsmatrix, Auswertung der Loesung bei T
        %wird in MonoVek geschrieben und diese zur Monodromiematrix
        %zusammengesetzt
        for k=1:Nz
            sol = ode45(@(psi,x)MathieuDGL(psi,x,D,nu_02,nu_C2),[t0,T],Diagonal(:,k));
            MonoVek = deval(sol,T);
            Monodromie(:,k) = MonoVek;
        end
    
        %charakteristische Multiplikatoren (Eigenwerte der Monodromiematrix)
        p = [1 -(Monodromie(1,1)+Monodromie(2,2)) (Monodromie(1,1)*Monodromie(2,2))-Monodromie(1,2)*Monodromie(2,1)];
        e = roots(p).';
        CharMult(l,:) = e;

        %Betrag der Eigenwerte zur Stabilitaetsanalyse
        eBetr = abs(e);
        
        Real = 1/T * log(e);
        

        %charakteristische Exponenten
        if nu_C2==nu_02
            CharEx(o,:) = [nu_02,nu_C2,Real(1,1),Real(1,2)];
            o=o+1;
        end
        

        % 1=stabil , -1=instabil Anmerkung nicht vergessen (zu =1)
        % sortieren der Parameterkombinationen in Array mit stabilen und
        % mit instabilen Kombinationen von nu_02 und nu_C2 anhand des
        % charakteristischen Mutiplikators
        b=0;

        if eBetr(1)<1 && eBetr(2)<1
            b=1;
            plotwertstabil(l,:) = [nu_02,nu_C2,b];
        elseif eBetr(1)>1 || eBetr(2)>1
            b=-1;
            plotwertinstabil(l,:) = [nu_02,nu_C2,b];
        end
    
        l=l+1;  
    end
end


%Loeschen der reinen 0-Zeilen
plotwertstabil = plotwertstabil(any(plotwertstabil,2),:);
plotwertinstabil = plotwertinstabil(any(plotwertinstabil,2),:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Grafische Darstellung: Struttsche Karte  (Plotten aller Wertepaare
%nu_02/nu_C2, die zu stabiler Loesung fuehren) sowie Realteile der
%charakteristischen Exponenten

figure('Position', [10 10 600 600])
tiledlayout(9,1);

nexttile([7 1])
% Struttsche Karte
scatter(plotwertstabil(:,1),plotwertstabil(:,2),1)
ylabel('$\rm{Parameter} \; \nu_C^2 \;\;\rm{[-]}$','interpreter','latex','Position', [-0.5 nu_C2/2],'FontSize', 10);
%xlabel('$\rm{Parameter} \; \nu_0^2 \;\;\rm{[-]}$','interpreter','latex','FontSize', 10);


% Realteile char. Exponenten
nexttile([2 1])
xachse = unt0:SW:ob0;
plot(xachse,CharEx(:,3:4))
xlabel('$\rm{Parameter} \; \nu_0^2 \;\;\rm{[-]}$','interpreter','latex','FontSize', 10);
ylabel('$Re(s_R) \;\; \rm{[-]}$','interpreter','latex','Position', [-0.5 -D],'FontSize', 10);


