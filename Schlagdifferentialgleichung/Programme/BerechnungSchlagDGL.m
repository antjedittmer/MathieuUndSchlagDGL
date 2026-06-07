%clear command window, clear workspace, close all figures
clc ; clear variables ; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialisierung und Hilfsvariablen

% true = 1;         % fuer Detektion von mu_param-Wert bei Trennung der Realtele
SW = 0.1;        % Schrittweite der Berechnung, Genauigkeit
idx = 1;            % Zaehlvariable
t0 = 0.0;         % Anfangszeitpunkt t0
T = 2*pi;         % Periodendauer
plotAll = 0;       % Diagnostic plots

excelDir = 'excelDir';
if ~isfolder(excelDir)
    mkdir(excelDir);
end

figDir = 'figDir';
if ~isfolder(figDir)
    mkdir(figDir);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Auswahl, welcher Rotor berechnet werden soll
%(bitte entkommentieren)

%Auswahl=1;Blatt=3; %3-Blatt-Rotor, see-saw
%Auswahl=2;Blatt=3; %3-Blatt-Rotor, voll gelenkig
%Auswahl=3;Blatt=4; %4-Blatt-Rotor, voll gelenkig
%Auswahl=4;Blatt=5; %5-Blatt-Rotor, voll gelenkig
%Auswahl=5;Blatt=3; %3-Blatt-Rotor, gelenk-/lagerlos
%Auswahl=6;Blatt=4; %4-Blatt-Rotor, gelenk-/lagerlos
%Auswahl = 7; Blatt = 1; %Einzelblattkoordinaten im rotierenden System

AuswahlInfo = {
    1, 3, '3-Blatt-Rotor, see-saw';
    2, 3, '3-Blatt-Rotor, voll gelenkig';
    3, 4, '4-Blatt-Rotor, voll gelenkig';
    4, 5, '5-Blatt-Rotor, voll gelenkig';
    5, 3, '3-Blatt-Rotor, gelenk-/lagerlos';
    6, 4, '4-Blatt-Rotor, gelenk-/lagerlos';
    7, 1, 'Einzelblattkoordinaten im rotierenden System'
    };


% if exist('Auswahl','var') ~= 1
%     Auswahl = 100;    % fuer Fehlermeldung, wenn keine Auswahl getroffen wurde
% end

%Exakte Berechung mittels Floquet oder Naeherung durch konstante
%Koeffizienten?

konstant=0;        %Exakt
%konstant=1;        %Naeherung

% if Auswahl == 100
%     f = warndlg('Bitte Rotorvariante waehlen!','Fehler');
%     return;
% end

for Auswahl = 1:3 %1: length(AuswahlInfo)

    % Reset counter and buffers for each Auswahl
    idx = 1;
    n = 1;
    cntN = 1;
    bufferDiffReNeg0 = 0;

    Blatt   = AuswahlInfo{Auswahl, 2};
    fprintf('Running: Auswahl %d — %s\n', Auswahl, AuswahlInfo{Auswahl, 3});

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    AnzGl = Blatt*2;      %Anzahl der Gleichungen
    b1 = AnzGl; %Blatt + 1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Laden der Datei mit Parametern des Rotors und der Berechnung

    Parameter = readtable('Parameter.xlsx','Range','C4:I29');
    Par = table2array(Parameter);

    rho = Par(26,1);

    ebeta = Par(8,Auswahl);
    gamma = Par(13,Auswahl);
    d2 = Par(17,Auswahl);
    d3 = Par(18,Auswahl);
    d4 = Par(19,Auswahl);
    nu0 = Par(20,Auswahl);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Grenzen fuer mu_param
    MuMin = 0;
    MuMax = 10;

    %%%
    ebetaStr = strrep(sprintf('_ebeta%2.3f', ebeta),'.','dot');
    muStr = sprintf('_mu%d_mu%d',round(MuMin), round(MuMax));
    auswahlStr = sprintf('_Auswahl%d',Auswahl);
    filename = ['Workspace', ebetaStr, muStr,auswahlStr,'.mat'];

    if exist(filename,'file') == 2

        load(filename, "damp","freq","MuMin","SW","MuMax", ...
            "CharExRe","CharExIm","nu0","Blatt","tableCharEx","CharExRe1Cor","CharExRe2Cor","CharExIm1Cor",'CharExIm2Cor');
    end

    if exist('CharExIm2Cor','var') ~= 1

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Vorbereitung der in den Schleifen zu fuellenden Arrays zwecks Programmbeschleunigung
        nMu = length(MuMin:SW:MuMax);
        Diagonal = diag(ones(AnzGl,1));
        Monodromie = zeros(AnzGl);
        CharMult = zeros(nMu,AnzGl+1);
        CharExRe = zeros(nMu,AnzGl);
        CharExIm = zeros(nMu,AnzGl);
        CharExImRaw = zeros(nMu,AnzGl);
        CharEx = zeros(nMu,AnzGl*6);

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
        buffer.Pos = zeros(Blatt,1);
        bufferDiffReNeg0 = 0; %
        cntN = 1;

        nuCSwitchVec = [4.4,5.7,8,9,11] - 0.15;
        %nuCSwitchValVec = [1.5,2.5,4.5,5,5];
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
            charMult = (eig(Monodromie))';
            [~,idxSort] = sort(real(charMult));
            charMultSort1 = charMult(idxSort);

            %Re = real(charMult);
            absMult1Greater = 0; %abs(real(charMultSort1(1))) - abs(real(charMultSort1(b1)))> 0.1;

            charMultSort = charMultSort1;
            prevAbsMult1Greater = absMult1Greater;

            CharMult(idx,:) = [charMultSort,mu_param];

            % charakteristische Exponenten
            if mu_param == 0
                Im0 = 1/T * angle(charMultSort);
            end

            % Berechne Real-und Imaginaerteile der Exponenten
            Eig.Real1 = 1/T * log(abs(charMultSort1));
            Eig.Real = 1/T * log(abs(charMultSort));
            Eig.Imag = 1/T * atan(imag(charMultSort)./real(charMultSort));

            % Korrigiere Imaginaerteil fuer kontinuierlichen
            % Verlauf
            % [Eig,buffer] = correctImagValues(Eig,buffer);

            for idxC = 1: length(Eig.Real)/2
                idxVec = 2*idxC-1 : 2*idxC;
                EigTemp.Real = Eig.Real(idxVec);
                EigTemp.Imag = Eig.Imag(idxVec);
                bufferTemp.Pos = buffer.Pos(idxC);

                [EigTemp,bufferTemp] = correctImagValues(EigTemp,bufferTemp);
                Eig.Real(idxVec) = EigTemp.Real;
                Eig.Imag(idxVec) = EigTemp.Imag;
                Eig.ImagSort(idxVec) = EigTemp.ImagSort;
                Eig.ImagCorrected(idxC) = EigTemp.ImagCorrected;
                Eig.ImagCorrectedNeg(idxC) = EigTemp.ImagCorrectedNeg;
                buffer.Pos(idxC) = bufferTemp.Pos;
            end

            %% Additionsterm n*2*pi/T
            absDiffCharMult = abs(CharMult(idx,1)) > abs(CharMult(idx,b1));
            if idx >1 && abs(prevAbsDiffCharMult - absDiffCharMult) > 0.1
                n =  n + 0.5;
                cntN = cntN+1;
            end
            prevAbsDiffCharMult = absDiffCharMult;

            %bufferDiffReNeg0 = diffReNeg0;

            nAdd = n*2*pi/T;
            nAddVec =  [-nAdd*ones(1,Blatt), nAdd*ones(1,Blatt)];
            ImagEigSortN = [Eig.ImagCorrectedNeg, Eig.ImagCorrected] + nAddVec;
            CharEx(idx,:) = [Eig.Real, Eig.Imag, Eig.ImagSort, ImagEigSortN, ...
                real(charMultSort), imag(charMultSort)];

            Re = Eig.Real;
            Im = Eig.Imag +  nAddVec; %ImagEigSortN; % sort(angle(eP(sortIdx)),'descend');

            CharExRe(idx,:) = Re;
            CharExIm(idx,:) = Im;
            CharExImRaw(idx,:) = Eig.Imag; %#ok<SAGROW>

            idx = idx+1;
        end

        %% Table fuer Excel Export
        lEig = length(Eig.Real);
        cellRealCharacteristicExponent = regexp(sprintf('RealCharExp%02d,', 1: length(Eig.Real)),',','split');
        cellImagCharacteristicExponent = regexp(sprintf('ImagCharExp%02d,', 1: length(Eig.Imag)),',','split');
        cellImagSorted = regexp(sprintf('ImagSort%02d,', 1: lEig),',','split');
        cellImagSortedN = regexp(sprintf('ImagSortN%02d,', 1: lEig),',','split');
        cellRealCharacteristicMultiplier = regexp(sprintf('RealCharMult%02d,', 1: lEig),',','split');
        cellImagCharacteristicMultiplier = regexp(sprintf('ImagCharMult%02d,', 1: lEig),',','split');

        % VarNames = {'muChar','RealCharExp1','RealCharExp2','ImagCharExp1','ImagCharExp2','minImagSort','maxImagSort',...
        %     'ImagSortN1','ImagSortN2','RealCharMult1','ImagCharMult1','RealCharMult2', 'ImagCharMult2'};
        VarNames = ['muChar',cellRealCharacteristicExponent(1:lEig),...
            cellImagCharacteristicExponent(1:lEig),...
            cellImagSorted(1:lEig), cellImagSortedN(1:lEig),...
            cellRealCharacteristicMultiplier(1:lEig),...
            cellImagCharacteristicMultiplier(1:lEig)];

        % Table mit allen Variablen (acuh sortierte Imaginaerwerte, die nur zum
        % Debuggen benutzt werden
        tableCharEx = array2table([CharMult(:,end),CharEx],"VariableNames",VarNames);

        % Variablen die ausgedruckt werden: Charakteristische Exponenten und
        % Multiplikatoren
        PrintNames = VarNames(contains(VarNames,'Char')); % Variablennamen
        tableCharPrint = tableCharEx(:,PrintNames ); % Auswahl der Tablespalten

        % 'Geglaettete' Realanteile
        CharExRe1 = CharExRe(:,1);
        CharExRe2 = CharExRe(:,Blatt + 1);

        CharExRePos = max(CharExRe1,CharExRe2); % Positiver Eigenwert: Glatter Verlauf
        CharExRe1_NegIdx =  abs(CharExRePos - CharExRe1) > eps; % Index: Negativer Wert 1. Realteil
        offset = CharExRe1(1); % Offset
        CharExReNeg =  - (CharExRePos - offset) + offset; % Negativer EW: Gespiegelter positiver EW

        CharExRe1Cor = CharExRe1; % Korrigierter 1. Eigenwert
        CharExRe1Cor(CharExRe1_NegIdx) = CharExReNeg(CharExRe1_NegIdx); % Ersetze negative Werte durch gespiegelte positive EW

        CharExRe2Cor = CharExRe2; % Korrigierter 1. Eigenwert
        CharExRe2Cor(~CharExRe1_NegIdx) = CharExReNeg(~CharExRe1_NegIdx); % Ersetze negative Werte durch gespiegelte positive EW

        tableCharPrint.RealCharExp1Corrected = CharExRe1Cor;
        tableCharPrint.RealCharExp2Corrected = CharExRe2Cor;

        % Raw imaginary parts (same ordering issues as real parts)
        CharExIm1 = CharExIm(:,1);
        CharExIm2 = CharExIm(:,Blatt+1);

        CharExIm1Cor = CharExIm1;
        CharExIm2Cor = CharExIm2;

        % At each mu step, assign imaginary parts to follow the larger/smaller
        % real part branch, consistent with how CharExRe1Cor/CharExRe2Cor were built
        for k = 1:length(mu_paramVec)
            if CharExRe1Cor(k) == CharExRe1(k)
                % Col1 is already the correct branch, no swap needed
                CharExIm1Cor(k) = CharExIm1(k);
                CharExIm2Cor(k) = -CharExIm1(k);
            else
                % Col1 was swapped in real part correction, so swap imaginary too
                CharExIm1Cor(k) = -CharExIm2(k);
                CharExIm2Cor(k) = CharExIm2(k);
            end
        end

        % Add to table for export
        tableCharPrint.ImagCharExp1Corrected = CharExIm1Cor;
        tableCharPrint.ImagCharExp2Corrected = CharExIm2Cor;

        excelfilename = strrep(filename,'.mat','.xlsx');
        excelfilename = strrep(excelfilename,'Workspace','CharactExponentenMultiplikatoren');

        save("Workspace.mat","damp","freq","MuMin","SW","MuMax", ...
            "CharExRe","CharExIm","nu0","Blatt","tableCharEx","CharExRe1Cor","CharExRe2Cor","CharExIm1Cor","CharExIm2Cor");
        save(filename, "damp","freq","MuMin","SW","MuMax", ...
            "CharExRe","CharExIm","nu0","Blatt","tableCharEx","CharExRe1Cor","CharExRe2Cor","CharExIm1Cor","CharExIm2Cor");

        excelfilename1 = fullfile(excelDir,excelfilename);
        writetable(tableCharPrint,excelfilename1);
    end


    %% Code for corrected and not corrected

    rotorDescription = AuswahlInfo{Auswahl, 3};

    figure(Auswahl);
    t = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
    title(t, sprintf('Auswahl: %d; Blatt: %d — %s', Auswahl, Blatt, rotorDescription));
    mu_vec = MuMin:SW:MuMax;
    cl = lines;

    % --- Subplot 1: Real parts (raw and corrected) ---
    ax2(1) = nexttile;
    plot(mu_vec, CharExRe1Cor, '-',  'Color', cl(1,:)); hold on;
    plot(mu_vec, CharExRe2Cor, '-',  'Color', cl(2,:));
    plot(mu_vec, CharExRe(:,1),  'k--', 'LineWidth', 1);
    plot(mu_vec, CharExRe(:,b1), 'k-.', 'LineWidth', 1);
    grid on;
    ylabel('Real part (-)');
    legend('Re(Exp1) cor', sprintf('Re(Exp%d) cor', Blatt), ...
        'Re(Exp1) raw',  sprintf('Re(Exp%d) raw', Blatt), ...
        'Location','northeastoutside');

    % --- Subplot 2: Imaginary parts (raw and corrected) ---
    ax2(2) = nexttile;
    plot(mu_vec, CharExIm1Cor, '-',  'Color', cl(1,:)); hold on;
    plot(mu_vec, CharExIm2Cor, '-',  'Color', cl(2,:));
    plot(mu_vec, CharExIm(:,1),  'k--', 'LineWidth', 1.0);
    plot(mu_vec, CharExIm(:,b1), 'k-.', 'LineWidth', 1.0);
    grid on;
    xlabel('Advance ratio \mu (-)');
    ylabel('Imaginary part (-)');
    legend('Im(Exp1) cor', sprintf('Im(Exp%d) cor', Blatt), ...
        'Im(Exp1) raw',  sprintf('Im(Exp%d) raw', Blatt), ...
        'Location','northeastoutside');

    linkaxes(ax2, 'x');

end

if plotAll == 0
    return;
end
%% Code von Matthieuscher DGL
figure(100);
ax1(1) = subplot(2,1,1);
plot(MuMin:SW:MuMax,CharExRe(:,1),'*-', MuMin:SW:MuMax,CharExRe(:,b1),'o-',...
    MuMin:SW:MuMax,CharExIm(:,1),'*-',MuMin:SW:MuMax,CharExIm(:,b1),'o-',...
    MuMin:SW:MuMax,CharExImRaw(:,1),'*-',MuMin:SW:MuMax,CharExImRaw(:,b1),'o-');
legend('Re(Exp1)', sprintf('Re(Exp%d)',b1),'Im(Exp1)',sprintf('Im(Exp%d)',b1),...
    'ImRaw(Exp1)',sprintf('ImRaw(Exp%d)',b1),...
    'Location','SouthWest'); grid on;
title(sprintf('Auswahl: %d; Blatt: %d', Auswahl, Blatt));
ax1(2) = subplot(2,1,2);
plot(MuMin:SW:MuMax,CharExRe1Cor,'*-', MuMin:SW:MuMax,CharExRe2Cor,'o-',...
    MuMin:SW:MuMax, CharExIm1Cor,'*-',MuMin:SW:MuMax,CharExIm2Cor,'o-');
grid on;
linkaxes(ax1,'x')

figure(102)
lineCell = {'-','--','-.','-','--','-.'};

for idx = 1: Blatt
    if Blatt > 1
        subplot(Blatt,1,idx);
    end
    plot(MuMin:SW:MuMax,CharExRe(:,idx),lineCell{idx},...
        MuMin:SW:MuMax,CharExRe(:,idx+Blatt),lineCell{idx+Blatt});
    grid on;
    legend(cellRealCharacteristicExponent([idx,idx+Blatt]))
    if idx == 1

        title(sprintf('Auswahl: %d; Blatt: %d', Auswahl, Blatt));
    end
end
xlabel('\mu (-)')



% evSort = real(CharMult(:,1:2))';
% figure;
% subplot(2,1,1);
% plot(mu_paramVec,eigenvalues_sorted(1,:), mu_paramVec,evSort(1,:),'k--');
% subplot(2,1,1);
% CharMult

figure(3);
plot(MuMin:SW:(MuMax-SW),abs(diff(CharExRe1Cor)),'-*')
legend('abs(diff(CharExRe1Cor))')

%% Berechnung absoluter Wert charakteristische Multiplikatoren
absCharMult1 = abs(CharMult(:,1));
absCharMultb1 = abs(CharMult(:,b1));
diffAbsCharMult = abs(absCharMult1) > abs(absCharMultb1);

checkCharMult1Greater = sum(absCharMult1 <= absCharMultb1) == length(absCharMult1);

vecAbsDiffRealCharMult = abs(diff(abs(absCharMult1) - abs(absCharMultb1)> 0));

figure(4);
pos0 = get(0,'defaultFigurePosition');
set(gcf, "Position",[pos0(1:3),1.5*pos0(4)])
ax1(1) = subplot(3,1,1);
plot(mu_paramVec,absCharMult1, mu_paramVec,absCharMultb1);
legend('Betrag char. Mult. 1',sprintf('Betrag char. Mult. %d',b1),'Location','SouthWest');
ylabel('linear')
axis tight; grid on;

ax1(2) = subplot(3,1,2);
semilogy(mu_paramVec,absCharMult1, mu_paramVec,absCharMultb1);
legend('Betrag char. Mult. 1',sprintf('Betrag char. Mult. %d',b1),'Location','SouthWest');
ylabel('semilogy')
axis tight; grid on;

ax1(3) = subplot(3,1,3);
plot(mu_paramVec,diffAbsCharMult,mu_paramVec(1:end-1),diff(diffAbsCharMult))
legend(sprintf('x = (|char. Mult. 1| > |char. Mult. %d|)',b1),'diff(x)','Location','SouthWest');
axis tight; grid on;
xlabel('\mu (-)');

% ax1(3) = subplot(3,1,3)
% plot(mu_paramVec(1:end-1),abs(diff(abs(absCharMult1) - abs(absCharMultb1)> 0)))


idxMuVec = mu_paramVec < 3;

mu_paramVec_idx = mu_paramVec(idxMuVec);

% figure;
% subplot(2,1,1)
% plot(mu_paramVec_idx,real(CharMult(idxMuVec,1)),...
%     mu_paramVec_idx,real(CharMult(idxMuVec,b1)));




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

