%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Berechnung der Struttschen Karten in den Grenzen von nu_02 und nu_C2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

loadMat = 0;  % mat-file laden, wenn Ergebnisse mit gleichem D vorhanden

SW = 0.1; %stepwidth
unt0 = 0; SW;
untC = 0; SW;
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
DVec = 0.001; %[0.15; 0.001; 0.2];

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
    plotwertstabil = zeros(lenNu,3);

    %nuCSwitchVec = [0.1,1,1.5^2,2^2,2.5^2] - 0.1;
    nuCSwitchVec = [0.2,1,1.5^2,2^2,2.5^2] - 0.15;
    n = 0;
    cntN = 1;

    buffer.Pos = 0;
    bufferNeg = 0;
    vecSwitch = zeros(lenNu02,1);

    noF = 1; % Number der Figure

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if exist(fileName,'file') == 2 && loadMat == 1
        load(fileName,'Char*','plotwert*', 'nu*');

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

                for k = 1 : Nz
                    sol = ode45(@(psi,x)MathieuDGL(psi,x,D,nu_02,nu_C2),[t0,T],Diagonal(:,k),options);
                    MonoVek  = deval(sol,T);
                    Monodromie(:,k) = MonoVek;
                end

                % Charakteristische Multiplikatoren (Eigenwerte der Monodromiematrix)
                p = [1, -(Monodromie(1,1) + Monodromie(2,2)), ...
                    Monodromie(1,1)*Monodromie(2,2)-Monodromie(1,2)*Monodromie(2,1)];
                eP = roots(p);
                %eBetr = abs(eP);

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
        save(fileName,'Char*','plotwert*', 'nu*');
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

end

%% Unused Code
% if SW <= 0.01 && 0
%     idxVec = 1:25:length(CharEx);
%     fileName1 = strrep(fileName,'SW1dot0e-02','SW5dot0e-02');
%     load(fileName1,'CharEx')
%     SW1 = 0.05;
% else
% idxVec = 1:length(CharEx);
% SW1 = SW;
%
% fileName1 = strrep(fileName,'.mat','_Diag.mat');
% tmp = load(fileName1);
% CharEx = tmp.CharEx;
%
%  save(fileName,'Char*','plotwert*', 'nu*');

%end
