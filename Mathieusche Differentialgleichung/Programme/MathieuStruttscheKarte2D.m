%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Berechnung der Struttschen Karten in den Grenzen von nu_02 und nu_C2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;
loadMat = 0;  % mat-file laden, wenn Ergebnisse mit gleichem D vorhanden
SW = 0.1; %stepwidth
unt0 = 0;
untC = 0;
ob0 = 9;
obC = 9;
fDir = 'figureFolder1'; % Ordner Abbildungen
if ~isdir(fDir) %#ok<ISDIR>
    mkdir(fDir)
end
dDir = 'dataFolder'; % Ordner mat-files
if ~isdir(dDir) %#ok<ISDIR>
    mkdir(dDir)
end
excelDir = 'dataFolder1';
if ~isdir(excelDir) %#ok<ISDIR>
    mkdir(excelDir);
end
% Anzahl der Gleichungen des DGL-Systems
Nz = 2;
% Parameter
DVec = 0.15; %0.001; %[0.15; 0.001; 0.2]; 0.3; %
% strDVec = {sprintf('%2.2f',DVec(1)); sprintf('%2.3f',DVec(2)); ...
%     sprintf('%2.1f',DVec(3))};
t0 = 0.0;
T = 2*pi;
for dIdx = 1: length(DVec)
    D = DVec(dIdx);
    % Anfangsbedingungen
    Diagonal = diag(ones(Nz,1));
    % Matlab mat-file Name
    matName = [strrep(sprintf('STRUTTscheKarte_D%2.1e_SW%2.1e',D,SW),'.','dot'),...
        '.mat'];
    if unt0 == 0
        matName = strrep(matName,'.mat','_unt0.mat');
    end
    fileName = fullfile(dDir,matName);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Vorbereitung der in den Schleifen zu fuellenden Arrays zwecks
    % Programmbeschleunigung
    lenNu02 = length(unt0:SW:ob0); % Anzahl Nu_02 Werte
    lenNuC2 = length(untC:SW:obC); % Anzahl Nu_C2 Werte
    lenNu = lenNu02 * lenNuC2; % Anzahl Kombination Nu_C2 Werte
    lenNuDiag = min(lenNu02,lenNuC2); % Anzahl Werte nu_02 == nu_C2
    Monodromie = zeros(Nz);
    CharEx = zeros(lenNuDiag,Nz*4+2); % nu, nc, Real1,2, Imag12, Imag1,2, Pole
    nAddVector = nan(lenNuDiag,1);
    plotwertstabil = zeros(lenNu,3);
    %nuCSwitchVec = [0.1,1,1.5^2,2^2,2.5^2] - 0.1;
    nuCSwitchVec = [0.2,1,1.5^2,2^2,2.5^2] - 0.15;
    n = 0;
    cntN = 1;
    buffer.Pos = 0;
    buffer.Neg = 0; % Added this line, as it was missing in the original code, but used in the function
    vecSwitch = zeros(lenNu02,1);
    noF = 1; % Number der Figure
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if exist(fileName,'file') == 2 && loadMat == 1
        load(fileName,'Char*','plotwert*', 'nu*','nAddVector');
    else
        % Zaehlvariablen
        lidx = 1; % Index charakteristische Multiplikatoren
        oidx = 1; % Index charakteristische Exponenten
        for nu_02 = unt0:SW:ob0
            for nu_C2  = untC:SW:obC
                %Schleife intergiert Nz Mal numerisch nacheinaner fuer die
                %Spaltenvektren der Einheitsmatrix, Auswertung der Loesung bei T
                %wird in MonoVek geschrieben und diese zur Monodromiematrix
                %zusammengesetzt
                options = odeset('RelTol',1e-10,'AbsTol',1e-12);
                
                % DGL System muss als separate function 'MathieuDGL(psi,x,D,nu_02,nu_C2)' existieren
                for k = 1 : Nz        
                        sol = ode45(@(psi,x)MathieuDGL(psi,x,D,nu_02,nu_C2),[t0,T],Diagonal(:,k),options);
                        MonoVek  = deval(sol,T);
                        Monodromie(:,k) = MonoVek;
                end
                
                % Characteristic Multipliers (Eigenwerte der Monodromiematrix)
                eP = eig(Monodromie);
                %% Charakteristische Exponenten
                if nu_C2 == nu_02
                    % Berechne Real-und Imaginaerteile der Exponenten
                    Eig.Real = 1/T * log(abs(eP));
                    Eig.Imag = 1/T * atan(imag(eP)./real(eP));
                    % Korrigiere Imaginaerteil fuer kontinuierlichen
                    % Verlauf
                    [Eig,buffer] = correctImagValues(Eig,buffer);
                    % % Additionsterm n*2*pi/T
                    if cntN <= length(nuCSwitchVec) && ... % Sicherheitscheck, damit n nicht groesser wird als die Anzahl der 'Switchstellen'
                            nu_C2 > nuCSwitchVec(cntN) && ... % n nur bei erreichen der Switchstellen umstellen
                            abs(Eig.Real(1) - Eig.Real(2))> eps % die Imaginaerteile muessen gleich sein
                        n =  n + 0.5;
                        cntN = cntN+1;
                    end
                    nAdd = n*2*pi/T; %
                    ImagEigSortN = [Eig.ImagCorrectedNeg, Eig.ImagCorrected] + [-nAdd,nAdd];
                    CharEx(oidx,:) = [nu_02, nu_C2, Eig.Real', min(Eig.ImagSort), max(Eig.ImagSort), ImagEigSortN, eP'];
                    nAddVector(oidx) = nAdd;
                    oidx = oidx + 1;
                end
                %% 1=stabile Kombinationen von nu_02 und nu_C2 basierend auf den
                % charakteristic Mutiplikators
                if  max(abs(eP)) < 1 %eBetr(1)<1 && eBetr(2)<1
                    b = 1;
                    plotwertstabil(lidx,:) = [nu_02,nu_C2,b];
                end
                lidx = lidx + 1;
            end
        end
        % Loeschen der reinen 0-Zeilen
        plotwertstabil =  plotwertstabil(any(plotwertstabil,2),:);
        save(fileName,'Char*','plotwert*', 'nu*','nAddVector');
    end
    %%
    try
        %[nu_02, nu_C2, Eig.Real', min(ImagEigSort), max(ImagEigSort), ImagEigSortN, eP'];
        CharExTable = array2table(CharEx(:,[1:4,7,8]),'VariableNames',...
            {'nu02','nu_C2','Eig.Real1','Eig.Real2', 'ImagEig1', 'ImagEig2'});
        CharExTable.RealCharExp1 = real(CharEx(:,9));
        CharExTable.ImagCharExp1 = imag(CharEx(:,9));
        CharExTable.RealCharExp2 = real(CharEx(:,10));
        CharExTable.ImagCharExp2 = imag(CharEx(:,10));
        % CharAllTable = array2table(CharEx(:,[1:4,7,8,9,10]),'VariableNames',...
        %    {'nu02','nu_C2','Eig.Real1','Eig.Real2', 'ImagEig1', 'ImagEig2','CharExp1','CharExp2'});
        excelfilename = strrep(fileName,'.mat','CharExAll.xlsx');
        excelfilename1 = strrep(excelfilename,dDir,excelDir);
        writetable(CharExTable,excelfilename1)
        % plotwertstabilTable = array2table(plotwertstabil,'VariableNames', {'nu02', 'nu_C2','b'});
        % excelfilename1 = strrep(fileName,'.mat','plotwertstabil.xlsx');
        % excelfilename1 = strrep(excelfilename1,'dataFolder','dataPlotWerteFolder');
        % if ~isdir('dataPlotWerteFolder') %#ok<ISDIR>
        %     mkdir('dataPlotWerteFolder')
        % end
        %
        % writetable(plotwertstabilTable, excelfilename1);
    catch
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Grafische Darstellung: Struttsche Karte  (Plotten aller Wertepaare
    %nu_02/nu_C2, die zu stabiler Loesung fuehren) sowie Realteile der
    % charakteristischen Exponenten
    cl = lines;
    fs = 10.5; %Fontsize
    hf = figure(dIdx);
    hf.Position = [10 10 600 600];
    try tiledlayout(11,1); catch, end% tiledlayout nicht unter R2007 verfuegbar
    try h(1) = nexttile([7 1]); catch, h(1) = subplot(11,1,[1,7]);  end %
    % Struttsche Karte
    scatter(plotwertstabil(:,1),plotwertstabil(:,2),5,'filled')
    ylabel('$\rm{Parameter} \; \nu_C^2 \;\;\rm{[-]}$','interpreter','latex','Position', [-0.5 obC/2],'FontSize', fs+2);
    grid on;
    % Realteile char. Exponenten
    try h(2) = nexttile([2 1]); catch, h(2) = subplot(11,1,[8,9]); end
    xachse = CharEx(:,1); %unt0:SW:ob0;
    plot(xachse,CharEx(:,3:4))
    ylabel('$Re(s_R) \;\; \rm{[-]}$','interpreter','latex','FontSize', fs); %'Position', [-0.5 -D]
    grid on;
    try h(3) = nexttile([2 1]); catch, h(3) = subplot(11,1,[10,11]); end %nexttile([2 1]);
    plot(xachse,CharEx(:,7:8)); hold on;
    % plot(xachse,CharEx(:,5),'--', 'color',cl(1,:)) % Zum Debuggen
    % plot(xachse,CharEx(:,6),'--', 'color',cl(2,:))
    grid on;
    ylim(gca,[-0.25,3]);
    grid on;
    xlabel('$\rm{Parameter} \; \nu_0^2 \;\;\rm{[-]}$','interpreter','latex','FontSize', fs+2);
    ylabel('$Im(s_R) \;\; \rm{[-]}$','interpreter','latex','FontSize', fs); %'Position', [-0.5 -D]
    linkaxes(h,'x')
    for idxH = 1: length(h)
        set(h(idxH),'TickLabelInterpreter','Latex','FontSize',fs)
    end
    pngname = fullfile(fDir,strrep(matName,'.mat',''));
    print(pngname, '-dpng')
    %
    hf1 = figure(dIdx+100);
    try tiledlayout(2,1); catch, end% tiledlayout nicht unter R2007 verfuegbar
    try h(1) = nexttile; catch, h(1) = subplot(2,1,1);  end %
    plot(xachse,CharEx(:,3:4)); grid on;
    title({['$\ddot{\phi} + 2D \dot{\phi} + (\nu^2_0 + \nu^2_C \cos(\psi))\phi = 0$, D = ', num2str(D,2), ', $\nu_C = \nu_0 $']...
        '$\Re(s_{Ri}) = 1/2\pi\ln(|\mu_{Ri}|)$, $\mu_{Ri}$: Eigenvalue monodromy matrix'},'interpreter','latex','FontSize', fs+2) % \text{atan}
    ylabel('$\Re(s_R) \;\; \rm{[-]}$','interpreter','latex','FontSize', fs+2);
    try h(1) = nexttile; catch, h(1) = subplot(2,1,2);  end %
    plot(xachse,CharEx(:,7:8)); hold on;  grid on;
    try
        %plot(xachse,nAddVector,'--','Color',0.5*ones(3,1))
        idxMVec = [0,diff(nAddVector)]>0;
        idxM = xachse(idxMVec);
        for idxPl = 1: length(idxM)
            plot(idxM(idxPl)*ones(2,1), get(gca,'ylim'),'Color',0.5*ones(3,1),'LineWidth',1);
        end
        for idxPl = 1: length(idxM)
            tmp = nAddVector(idxMVec);
            text( idxM(idxPl) + 0.1, max(get(gca,'ylim'))-1, num2str(tmp(idxPl),2),'interpreter','latex','FontSize', fs+1);
        end
        nAddVector(idxMVec)
        legend('$s_{R1}$','$s_{R2}$','increase in m', 'interpreter','latex','Location','SouthEast','FontSize', fs+1)
        title('$\Im(s_{Ri}) = 1/2\pi (\arctan(\Im(\mu_{Ri})/\Re(\mu_{Ri})) + m$','interpreter','latex','FontSize', fs+2)
    catch
    end
    xlabel('$\rm{Parameter} \; \nu_C^2 \;\;\rm{[-]}$','interpreter','latex','FontSize', fs+2);
    ylabel('$\Im(s_R) \;\; \rm{[-]}$','interpreter','latex','FontSize', fs+2); %'Position', [-0.5 -D]
    pngname = fullfile(fDir,strrep(strrep(strrep(matName,'.mat',''),'STRUTTscheKarte','CharExp'),'_unt0',''));
    print(pngname, '-dpng')
end
%%
% This function was part of the original code snippet (appended at the end)
% and is required by the main script for calculating characteristic exponents.
function  [Eig, buffer] = correctImagValues(Eig, buffer)
% correctImagValues korrigiert die Imaginaerwerte fuer stetigen Verlauf
% Inputs
% - Eig: Struct mit Vektor Eig.Imag mit zwei Werten der Imaginaerteile
% - buffer struct: Puffer mit letzten Werten fuer Maximum und Minimum
% Outputs
% - Eig: Struct mit angehaengtem, korrigierten Imaginaerteilen
% - buffer struct: Puffer ueberschrieben mit neuen Werten fuer Maximum und Minimum
% Zwei Checks fuer das korrekte Format
if nargin~= 2
    error('Two inputs are expected: The current imaginary part of the eigenvalues and the buffer with the last values');
end
if max(size(Eig.Imag))~= 2 || min(size(Eig.Imag))~= 1
    error('The current imaginary parts of an eigenvalue pair is expected');
end
% Sortiere Imaginaerteil
Eig.ImagSort = sort(Eig.Imag);
% Imaginaeranteil kontinuierlich steigend oder fallend
tmp = Eig.ImagSort(2);
tmpNeg = Eig.ImagSort(1);

% Initialize buffer.Pos/buffer.Neg if they don't exist, though they should be
% initialized in the main loop to 0
if ~isfield(buffer, 'Pos'), buffer.Pos = 0; end
if ~isfield(buffer, 'Neg'), buffer.Neg = 0; end

if buffer.Pos <= tmp || (abs(tmp) < 10^-5) % Wert uebernehmen
    Eig.ImagCorrected = tmp; % steigender pos. Wert
    Eig.ImagCorrectedNeg = tmpNeg;  % fallender neg. Wert
    buffer.Pos = 0; % Removed based on original logic, these were only placeholders
    buffer.Neg = 0;
else  % 'korrigierten' Wert nehmen fuer kontinuierlichen Verlauf
    % Wird der 'else'-Zweig getriggert, ist der buffer
    % bei 90 Grad (1.57 rad). Der positive Imaginaerteil ist damit
    % die Summe aus 180 Grad und dem negativen Winkel, dessen
    % Wert von -90 Grad zu 0 Grad laeuft.
    Eig.ImagCorrected = 2*buffer.Pos  + tmpNeg; % korrigierter pos. Wert
    Eig.ImagCorrectedNeg  = 2*buffer.Neg + tmp; % korrigierter neg. Wert
end
% Buffer ueberschreiben mit aktuellem Wert
buffer.Pos = max(tmp,buffer.Pos); % Maximum speichern
buffer.Neg = min(tmpNeg,buffer.Neg); % Minimum speichern
% Anmerkung: hier kammen leider andere Werte raus als bei atanh
% Eig.ImagAngle = 1/T * angle(imag(eP)./real(eP));
end
