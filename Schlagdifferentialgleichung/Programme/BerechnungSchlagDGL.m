%clear command window, clear workspace, close all figures
clc ; clear variables ; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialisierung und Hilfsvariablen

% true = 1;         % fuer Detektion von mu_param-Wert bei Trennung der Realtele
SW = 0.01;        % Schrittweite der Berechnung, Genauigkeit
idx = 1;            % Zaehlvariable
t0 = 0.0;         % Anfangszeitpunkt t0
T = 2*pi;         % Periodendauer

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Auswahl, welcher Rotor berechnet werden soll
%(bitte entkommentieren)

%Auswahl=1;Blatt=3; %3-Blatt-Rotor, see-saw
%Auswahl=2;Blatt=3; %3-Blatt-Rotor, voll gelenkig
%Auswahl=3;Blatt=4; %4-Blatt-Rotor, voll gelenkig
%Auswahl=4;Blatt=5; %5-Blatt-Rotor, voll gelenkig
%Auswahl=5;Blatt=3; %3-Blatt-Rotor, gelenk-/lagerlos
%Auswahl=6;Blatt=4; %4-Blatt-Rotor, gelenk-/lagerlos
Auswahl = 7; Blatt=1; %Einzelblattkoordinaten im rotierenden System

if exist('Auswahl','var') ~= 1
    Auswahl = 100;    % fuer Fehlermeldung, wenn keine Auswahl getroffen wurde
end

%Exakte Berechung mittels Floquet oder Naeherung durch konstante
%Koeffizienten?

konstant=0;        %Exakt
%konstant=1;        %Naeherung

if Auswahl == 100
    f = warndlg('Bitte Rotorvariante waehlen!','Fehler');
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AnzGl = Blatt*2;      %Anzahl der Gleichungen

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Laden der Datei mit Parametern des Rotors und der Berechnung

Parameter = readtable('Parameter.xlsx','Range','C4:I29');
Par = table2array(Parameter);

rho = Par(26,1);

ebeta = 0.3; %Par(8,Auswahl);
gamma = Par(13,Auswahl);
d2 = Par(17,Auswahl);
d3 = Par(18,Auswahl);
d4 = Par(19,Auswahl);
nu0 = Par(20,Auswahl);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Grenzen fuer mu_param
MuMin = 0;
MuMax = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Vorbereitung der in den Schleifen zu fuellenden Arrays zwecks Programmbeschleunigung
nMu = length(MuMin:SW:MuMax);
Diagonal = diag(ones(AnzGl,1));
Monodromie = zeros(AnzGl);
CharMult = zeros(nMu,AnzGl+1);
CharExRe = zeros(nMu,AnzGl);
CharExIm = zeros(nMu,AnzGl);
CharEx = zeros(nMu,AnzGl*5 + 2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Einlesen der A-Matrix mit mu_param=0 fuer die exakte Loesung im Schwebeflug

[~,A] = SchlagDGL(0,0,gamma,d2,d3,d4,0,ebeta,nu0,Blatt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Berechnung der Eigenwerte von mu_min bis mu_max
%Berechnen der Eigenwerte im Schwebeflug
[freq,sortIdx] = sort(imag(eig(A)).','descend'); %#ok<UDIM> geringer Geschwindigkeitsverlust durch kleine Matrix
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


% Optionen fuer die Genauigkeit und Toleranz fuer den ODE-Solver
options1 = odeset('RelTol',1e-10,'AbsTol',1e-12);
options2 = odeset('RelTol',1e-10,'AbsTol',1e-12,'MaxStep', 1e-3);
options3 = odeset('RelTol', 1e-10, 'AbsTol', 1e-12, ...
    'Jacobian', @(psi, x) myJacobian(psi, x, gamma, mu_param, d2, d3, ebeta));

% 'Switch'-Stellen des Imaginaerteil
buffer.Pos = 0;
bufferDiffReNeg0 = 0; % 
cntN = 1;

nuCSwitchVec = [4.4,5.5,8,9,11] - 0.15;
nuCSwitchValVec = [1.5,2.5,4.5,5,5];
n = 1;

mu_paramVec = MuMin:SW:MuMax;
cond_A = nan(length(mu_paramVec),1); 

for mu_param = mu_paramVec

    if konstant == 1
        for k=1:AnzGl
            sol = ode45(@(psi,x)SchlagDGLkonstant(psi,x,gamma,d2,d3,d4,mu_param,ebeta,nu0,Blatt),[t0,T],Diagonal(:,k),options);
            MonoVek = deval(sol,T);
            Monodromie(:,k) = MonoVek;
        end
    elseif konstant == 0

        for k=1:AnzGl
            sol = ode45(@(psi, x) SchlagDGL(psi, x, gamma, d2, d3, d4, mu_param, ebeta, nu0, Blatt), ...
                [t0, T], Diagonal(:,k), options2);
            MonoVek = deval(sol,T);
            Monodromie(:,k) = MonoVek;
        end

        cond_A(idx) = cond(Monodromie);
    end

    % charakteristische Multiplikatoren (Eigenwerte der Monodromiematrix)
    eP = (eig(Monodromie))';
    ePSort = sort(eP); 
    CharMult(idx,:) = [ePSort,mu_param];

    % charakteristische Exponenten
    if mu_param == 0
        Im0 = 1/T * angle(eP(sortIdx));
    end

    % Berechne Real-und Imaginaerteile der Exponenten
    Eig.Real = 1/T * log(abs(ePSort));
    Eig.Imag = 1/T * atan(imag(ePSort)./real(ePSort));

    % Korrigiere Imaginaerteil fuer kontinuierlichen
    % Verlauf
    [Eig,buffer] = correctImagValues(Eig,buffer);

    % % Additionsterm n*2*pi/T

    diffReNeg0 = (Eig.Real(2) - Eig.Real(1)) >= 0 & abs(Eig.Real(2) - Eig.Real(1)) > 0.1;
    if cntN <= length(nuCSwitchVec) && ... % Sicherheitscheck, damit n nicht groesser wird als die Anzahl der 'Switchstellen'
            mu_param > nuCSwitchVec(cntN) && ... % n nur bei erreichen der Switchstellen umstellen
            abs(Eig.Imag(1) - Eig.Imag(2)) < eps && ...% die Imaginaerteile muessen gleich sein
            diffReNeg0 ~= bufferDiffReNeg0 % die Realteile 'kreuzen' sich
            % abs(Eig.Real(1) - Eig.Real(2)) < nuCSwitchValVec(cntN) % die Realteile 'kreuzen' sich
        n =  n + 0.5;
        cntN = cntN+1;
    end
    bufferDiffReNeg0 = diffReNeg0;

    nAdd = n*2*pi/T;
    ImagEigSortN = [Eig.ImagCorrectedNeg, Eig.ImagCorrected] + [-nAdd, nAdd];
    CharEx(idx,:) = [Eig.Real,  Eig.Imag, min(Eig.ImagSort), max(Eig.ImagSort), ImagEigSortN, ...
        real(ePSort(1)), imag(ePSort(1)), real(ePSort(2)), imag(ePSort(2))];

    Re = Eig.Real;
    Im = sort(Eig.Imag) + [-nAdd, nAdd]; %ImagEigSortN; % sort(angle(eP(sortIdx)),'descend');

    CharExRe(idx,:) = Re;
    CharExIm(idx,:) = Im;

    idx = idx+1;
end

%% Table fuer Excel Export

VarNames = {'muChar','RealCharExp1','RealCharExp2','ImagCharExp1','ImagCharExp2','minImagSort','maxImagSort',...
    'ImagSortN1','ImagSortN2','RealCharMult1','ImagCharMult1','RealCharMult2', 'ImagCharMult2'};

% Table mit allen Variablen (acuh sortierte Imaginaerwerte, die nur zum
% Debuggen benutzt werden
tableCharEx = array2table([CharMult(:,end),CharEx],"VariableNames",VarNames);

% Variablen die ausgedruckt werden: Charakteristische Exponenten und
% Multiplikatoren
PrintNames = VarNames(contains(VarNames,'Char')); % Variablennamen
tableCharPrint = tableCharEx(:,PrintNames ); % Auswahl der Tablespalten

% 'Geglaettete' Realanteile 
CharExRe1 = CharExRe(:,1);
CharExRe2 = CharExRe(:,2);

CharExRePos = max(CharExRe1,CharExRe2); % Positiver Eigenwert: Glatter Verlauf
CharExRe1_NegIdx =  abs(CharExRePos - CharExRe1) > eps; % Index: Negativer Wert 1. Realteil 
offset = CharExRe1(1); % Offset 
CharExReNeg =  - (CharExRePos - offset) + offset; % Negativer EW: Gespiegelter positiver EW

CharExRe1Cor = CharExRe1; % Korrigierter 1. Eigenwert
CharExRe1Cor(CharExRe1_NegIdx) = CharExReNeg(CharExRe1_NegIdx); % Ersetze negative Werte durch gespiegelte positive EW

CharExRe2Cor = CharExRe2; % Korrigierter 1. Eigenwert
CharExRe2Cor(~CharExRe1_NegIdx) = CharExReNeg(~CharExRe1_NegIdx); % Ersetze negative Werte durch gespiegelte positive EW


tableCharPrint.RealCharExp1Corrected = CharExRe1Cor;
tableCharPrint.ImagCharExp1Corrected = CharExIm(:,1); % Erster Wert 
tableCharPrint.RealCharExp2Corrected = CharExRe2Cor;
tableCharPrint.ImagCharExp2Corrected = - CharExIm(:,1);

excelDir = 'excelDir';
if ~isfolder(excelDir)
    mkdir(excelDir);
end

ebetaStr = strrep(sprintf('_ebeta%2.3f', ebeta),'.','dot');
muStr = sprintf('_mu%d_mu%d',round(MuMin), round(MuMax));
filename = ['Workspace', ebetaStr, muStr,'.mat'];

excelfilename = strrep(filename,'.mat','.xlsx');
excelfilename = strrep(excelfilename,'Workspace','CharactExponentenMultiplikatoren');

save("Workspace.mat","damp","freq","MuMin","SW","MuMax", ...
    "CharExRe","CharExIm","nu0","Blatt","tableCharEx");
save(filename, "damp","freq","MuMin","SW","MuMax", ...
    "CharExRe","CharExIm","nu0","Blatt","tableCharEx");

excelfilename1 = fullfile(excelDir,excelfilename);

writetable(tableCharPrint,excelfilename1);

%% Code von Matthieuscher DGL
% %[nu_02, nu_C2, Eig.Real', min(ImagEigSort), max(ImagEigSort), ImagEigSortN, eP'];
% CharExTable = array2table(CharEx(:,[1:4,7,8]),'VariableNames',...
%     {'nu02','nu_C2','Eig.Real1','Eig.Real2', 'ImagEig1', 'ImagEig2'});
% 
% CharExTable.RealCharExp1 = real(CharEx(:,9));
% CharExTable.ImagCharExp1 = imag(CharEx(:,9));
% CharExTable.RealCharExp2 = real(CharEx(:,10));
% CharExTable.ImagCharExp2 = imag(CharEx(:,10));
% 
% % CharAllTable = array2table(CharEx(:,[1:4,7,8,9,10]),'VariableNames',...
% %    {'nu02','nu_C2','Eig.Real1','Eig.Real2', 'ImagEig1', 'ImagEig2','CharExp1','CharExp2'});
% 
% 
% excelfilename = strrep(fileName,'.mat','CharExAll.xlsx');
% excelfilename1 = strrep(excelfilename,dDir,excelDir);
% 
% writetable(CharExTable,excelfilename1)


figure; plot(MuMin:SW:MuMax,CharExIm(:,1), MuMin:SW:MuMax,CharExIm(:,2));
hold on; plot(MuMin:SW:MuMax,CharExRe(:,1), MuMin:SW:MuMax, CharExRe(:,2));

% figure;
% plot(real(CharMult(:,1)),imag(CharMult(:,1)),'kx')
% hold on;
% plot(real(CharMult(:,2)),imag(CharMult(:,2)),'bx')
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Korrektur der Ergebnisse

% CharExIm = (1/T) * unwrap(CharExIm);

% for k=1:length(MuMin:SW:MuMax)
%    CharExIm(k,:) = CharExIm(k,:) + freqInt;
% end


%nur entkommentieren, wenn im gewählten Bereich für große mu_param Probleme
%auftreten: Schleife spiegelt oberen Ast der Realteile auf den unterern
%(Annahme der Symmetrie) und ersetzt fehlerhafte Imaginärteile durch
%Im(s_R(mu_param)) an der Stelle des gewählten mu_param

%{

%nu anpassen ab
mu_param = 1.2;
format shortG
nutemp = round(mu_param/SW);

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

