%clear command window, clear workspace, close all figures
clc ; clear variables ; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialisierung und Hilfsvariablen

Auswahl=100;    %fuer Fehlermeldung, wenn keine Auswahl getroffen wurde
true=1;         %fuer Detektion von mu-Wert bei Trennung der Realtele
SW = 0.01;      %Schrittweite der Berechnung, Genauigkeit
l=1;            %Zaehlvariable
t0 = 0.0;       %Anfangszeitpunkt t0
T = 2*pi;       %Periodendauer

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Auswahl, welcher Rotor berechnet werden soll
%(bitte entkommentieren)

Auswahl=1;Blatt=3; %3-Blatt-Rotor, see-saw
%Auswahl=2;Blatt=3; %3-Blatt-Rotor, voll gelenkig
%Auswahl=3;Blatt=4; %4-Blatt-Rotor, voll gelenkig
%Auswahl=4;Blatt=5; %5-Blatt-Rotor, voll gelenkig
%Auswahl=5;Blatt=3; %3-Blatt-Rotor, gelenk-/lagerlos
%Auswahl=6;Blatt=4; %4-Blatt-Rotor, gelenk-/lagerlos
%Auswahl=7;Blatt=1; %Einzelblattkoordinaten im rotierenden System

%Exakte Berechung mittels Floquet oder Naeherung durch konstante
%Koeffizienten?
 
konstant=0;        %Exakt
%konstant=1;        %Naeherung

if Auswahl == 100
    f = warndlg('Bitte Rotorvariante waehlen!','Fehler');
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AnzGl=Blatt*2;      %Anzahl der Gleichungen

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Laden der Datei mit Parametern des Rotors und der Berechnung

Parameter = readtable('Parameter.xlsx','Range','C4:I29');
Par = table2array(Parameter);

rho = Par(26,1);

ebeta=Par(8,Auswahl);
gamma=Par(13,Auswahl);
d2=Par(17,Auswahl);
d3=Par(18,Auswahl);
d4=Par(19,Auswahl);
nu0=Par(20,Auswahl);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Grenzen fuer mu
MuMin = 0;
MuMax = 4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Vorbereitung der in den Schleifen zu fuellenden Arrays zwecks Programmbeschleunigung

Diagonal=diag(ones(AnzGl,1));
Monodromie = zeros(AnzGl);
CharMult = zeros(length(MuMin:SW:MuMax),AnzGl+1);
CharExRe = zeros(length(MuMin:SW:MuMax),AnzGl);
CharExIm = zeros(length(MuMin:SW:MuMax),AnzGl);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Einlesen der A-Matrix mit mu=0 fuer die exakte Loesung im Schwebeflug

[~,A] = SchlagDGL(0,0,gamma,d2,d3,d4,0,ebeta,nu0,Blatt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Berechnung der Eigenwerte von mu_min bis mu_max
%Berechnen der Eigenwerte im Schwebeflug
[freq,sortIdx] = sort(imag(eig(A)).','descend');
damp = real(eig(A));

if freq(1)-fix(freq(1))<0.5
    freqInt = round(freq);
elseif freq(1)-fix(freq(1))>=0.5 && freq(1)-fix(freq(1))<1
    freqInt = -round(freq);
end

if Blatt == 5
    sortIdx2 = [1,2,3,6,7,4,5,8,9,10];
    freqInt = freqInt(sortIdx2);
end


%Optionen fuer die Genauigkeit und Toleranz fuer den ODE-Solver
options = odeset('RelTol',1e-10,'AbsTol',1e-12);

for m = MuMin:SW:MuMax
    mu = m ;
    if konstant == 1
        for k=1:AnzGl
            sol = ode45(@(psi,x)SchlagDGLkonstant(psi,x,gamma,d2,d3,d4,mu,ebeta,nu0,Blatt),[t0,T],Diagonal(:,k),options);
            MonoVek = deval(sol,T);
            Monodromie(:,k) = MonoVek;
        end
    elseif konstant == 0
        for k=1:AnzGl
            sol = ode45(@(psi,x)SchlagDGL(psi,x,gamma,d2,d3,d4,mu,ebeta,nu0,Blatt),[t0,T],Diagonal(:,k),options);
            MonoVek = deval(sol,T);
            Monodromie(:,k) = MonoVek;
        end
    end

    %charakteristische Multiplikatoren (Eigenwerte der Monodromiematrix)
    e = (eig(Monodromie)).';
    CharMult(l,:) = [e,mu];


    %charakteristische Exponenten
    if mu==0
        Im0 = 1/T * angle(e(sortIdx))
    end
    
        
    if freqInt==round(freq)
        Re = sort(1/T * log(abs(e)),'descend');
    else
        Re = sort(1/T * log(abs(e)),'ascend');
    end
    
    Im = sort(angle(e(sortIdx)),'descend');
        

    CharExRe(l,:) = Re;
    CharExIm(l,:) = Im;

    l=l+1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Korrektur der Ergebnisse

CharExIm = (1/T) * unwrap(CharExIm);

for k=1:length(MuMin:SW:MuMax)
    CharExIm(k,:) = CharExIm(k,:) + freqInt;
end


%nur entkommentieren, wenn im gewählten Bereich für große mu Probleme
%auftreten: Schleife spiegelt oberen Ast der Realteile auf den unterern
%(Annahme der Symmetrie) und ersetzt fehlerhafte Imaginärteile durch
%Im(s_R(mu)) an der Stelle des gewählten mu

%{

%nu anpassen ab
mu = 1.2;
format shortG
nutemp = round(mu/SW);

temp = CharExRe(1,1);
CharExRe=CharExRe-temp;
for k=nutemp:length(CharExRe)
    CharExRe(k,1:3) = -CharExRe(k,4:6);
end
CharExRe=CharExRe+temp;

for k=nutemp:length(CharExIm)
    CharExIm(k,:) = CharExIm(nutemp,:);
end
%}

save('Workspace.mat')