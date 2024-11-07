%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Berechnung der Struttschen Karten in den Grenzen von nu_02 und nu_C2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

loadMat = 0;  % mat-file laden, wenn Ergebnisse mit gleichem D vorhanden

SW = 0.01; %stepwidth
unt0 = 0; SW;
untC = 0; SW;
ob0 = 9;
obC = 9;

fDir = 'figureFolder'; % Ordner Abbildungen
if ~isdir(fDir) %#ok<ISDIR>
    mkdir(fDir)
end

dDir = 'dataFolder'; % Ordner mat-files
if ~isdir(dDir) %#ok<ISDIR>
    mkdir(dDir)
end

% Anzahl der Gleichungen des DGL-Systems
Nz = 2;

% Parameter
DVec = [0.15;0.001; 0.2]; 

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
    matName = strrep(matName,'.mat','_Diag.mat');

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

    nuCSwitchVec = [0.2,1,1.5^2,2^2,2.5^2] - 0.15;
    n = 0;
    cntN = 1;

    buffer = 0;
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

        %         % Charakteristische Exponenten
        % A = MathieuDGLlin(D,nu_02,nu_C2);
        % RealEig1(:,lidx) = real(eig(A)); %1/T * log(eP);
        % ImagEig1(:,lidx) = imag(eig(A)); %1/T * atan(imag(eP)./real(eP));

        for nu_02 = unt0:SW:ob0

            nu_C2 = nu_02;
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
            eBetr = abs(eP);

            %% Charakteristische Exponenten

            RealEig = 1/T * log(eBetr);
            ImagEig = 1/T * atan(imag(eP)./real(eP));
            [~,idx] = sort(ImagEig);
            ImagEigSort = ImagEig(idx);

            if cntN <= length(nuCSwitchVec) && ...
                    nu_C2 > nuCSwitchVec(cntN) && ...
                    abs(RealEig(1) - RealEig(2))> eps
                n =  n + 0.5;
                cntN = cntN+1;
            end

            % Imaginaeranteil kontinuierlich steigend oder fallend
            tmp = ImagEigSort(2);
            tmpNeg = ImagEigSort(1);

            if buffer <= tmp || (abs(tmp) < 10^-5) % Wert uebernehmen
                ImagEigCorrected = tmp;
                ImagEigCorrectedNeg = tmpNeg;
                buffer = 0;
                bufferNeg = 0;
            else % 'anderen' Verlauf nehmen fuer kontinuierlich
                ImagEigCorrected = 2*buffer + tmpNeg; % steigend oder fallend
                ImagEigCorrectedNeg = 2*bufferNeg + tmp;
                vecSwitch(oidx)= 1;

            end
            buffer = max(tmp,buffer); % Maximum speichern
            bufferNeg = min(tmpNeg,bufferNeg); % Minimum speichern

            % Additionsterm  n*2*pi/T
            nAdd = n*2*pi/T;
            ImagEigSortN = [ImagEigCorrectedNeg, ImagEigCorrected] + [-nAdd,nAdd];

            CharEx(oidx,:) = [nu_02, nu_C2, RealEig', min(ImagEigSort), max(ImagEigSort), ImagEigSortN, eP'];
            oidx = oidx + 1;

            %fprintf('LinEig: %2.2f,%2.2f, Eig: %2.2f,%2.2f\n',RealEig1(1),RealEig1(2),RealEig(1),RealEig(2))


            [~,idx] = sort(ImagEig);
            ImagEigSort = ImagEig(idx);

            if cntN <= length(nuCSwitchVec) && ...
                    nu_C2 > nuCSwitchVec(cntN) && ...
                    abs(RealEig(1) - RealEig(2))> eps
                n =  n + 0.5;
                cntN = cntN+1;
            end

            % Imaginaeranteil kontinuierlich steigend oder fallend
            tmp = ImagEigSort(2);
            tmpNeg = ImagEigSort(1);

            if buffer <= tmp || (abs(tmp) < 10^-5) % Wert uebernehmen
                ImagEigCorrected = tmp;
                ImagEigCorrectedNeg = tmpNeg;
                buffer = 0;
                bufferNeg = 0;
            else % 'anderen' Verlauf nehmen fuer kontinuierlich
                ImagEigCorrected = 2*buffer + tmpNeg; % steigend oder fallend
                ImagEigCorrectedNeg = 2*bufferNeg + tmp;
                vecSwitch(oidx)= 1;

            end
            buffer = max(tmp,buffer); % Maximum speichern
            bufferNeg = min(tmpNeg,bufferNeg); % Minimum speichern

            % Additionsterm  n*2*pi/T
            nAdd = n*2*pi/T;
            ImagEigSortN = [ImagEigCorrectedNeg, ImagEigCorrected] + [-nAdd,nAdd];

            CharEx(lidx,:) = [nu_02, nu_C2, RealEig', min(ImagEigSort), max(ImagEigSort), ImagEigSortN, eP'];  %
            %CharEx1(lidx,:) = [nu_02, nu_C2, RealEig', RealEig1', min(ImagEigSort), max(ImagEigSort), ImagEigSortN];

            % 1=stabile Kombinationen von nu_02 und nu_C2 basierend auf den
            % charakteristic Mutiplikators

            if eBetr(1)<1 && eBetr(2)<1
                b = 1;
                plotwertstabil(lidx,:) = [nu_02,nu_C2,b];
            end
            lidx = lidx + 1;


        end

        % Loeschen der reinen 0-Zeilen
        plotwertstabil =  plotwertstabil(any(plotwertstabil,2),:);
        save(fileName,'Char*','plotwert*', 'nu*');
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Grafische Darstellung: Struttsche Karte  (Plotten aller Wertepaare
    %nu_02/nu_C2, die zu stabiler Loesung fuehren) sowie Realteile der
    % charakteristischen Exponenten
    cl = lines;
    fs = 10.5; %Fontsize
    hf = figure(dIdx);
    hf.Position = [10 10 600 600];
    tiledlayout(11,1);


    h(1) = nexttile([7 1]);
    % Struttsche Karte
    scatter(plotwertstabil(:,1),plotwertstabil(:,2),5,'filled')
    ylabel('$\rm{Parameter} \; \nu_C^2 \;\;\rm{[-]}$','interpreter','latex','Position', [-0.5 obC/2],'FontSize', fs+2);
    grid on;

    % if SW <= 0.01
    %     idxVec = 1:25:length(CharEx);
    %     fileName1 = strrep(fileName,'SW1dot0e-02','SW5dot0e-02');
    %     load(fileName1,'CharEx')
    %     SW1 = 0.05;
    % else
    idxVec = 1:length(CharEx);
    SW1 = SW;
    % end

    % Realteile char. Exponenten
    h(2) = nexttile([2 1]);
    xachse = unt0:SW1:ob0;
    plot(xachse,CharEx(:,3:4))

    ylabel('$Re(s_R) \;\; \rm{[-]}$','interpreter','latex','FontSize', fs); %'Position', [-0.5 -D]
    grid on;

    h(3) = nexttile([2 1]);
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
