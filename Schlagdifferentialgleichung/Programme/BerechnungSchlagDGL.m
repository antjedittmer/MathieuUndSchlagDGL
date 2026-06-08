clc; clear variables; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialisierung

SW     = 0.1;
t0     = 0.0;
T      = 2*pi;
plotAll = 0;
loadMat = 0;

excelDir = 'excelDir';
if ~isfolder(excelDir); mkdir(excelDir); end

figDir = 'figDir';
if ~isfolder(figDir); mkdir(figDir); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Rotor-Auswahl

AuswahlInfo = {
    1, 3, '3-Blatt-Rotor, see-saw';
    2, 3, '3-Blatt-Rotor, voll gelenkig';
    3, 4, '4-Blatt-Rotor, voll gelenkig';
    4, 5, '5-Blatt-Rotor, voll gelenkig';
    5, 3, '3-Blatt-Rotor, gelenk-/lagerlos';
    6, 4, '4-Blatt-Rotor, gelenk-/lagerlos';
    7, 1, 'Einzelblattkoordinaten im rotierenden System'
    };


for Auswahl = 7

    %% Reset
    clearvars -except AuswahlInfo SW t0 T plotAll excelDir figDir Auswahl loadMat

    Blatt  = AuswahlInfo{Auswahl, 2};
    fprintf('Running: Auswahl %d — %s\n', Auswahl, AuswahlInfo{Auswahl, 3});

    AnzGl  = Blatt * 2;
    b1     = AnzGl;
    nPairs = Blatt;
    konstant = 0;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Parameter laden

    Parameter = readtable('Parameter.xlsx', 'Range', 'C4:I29');
    Par       = table2array(Parameter);

    rho   = Par(26,1);
    ebeta = Par(8,  Auswahl);
    gamma = Par(13, Auswahl);
    d2    = Par(17, Auswahl);
    d3    = Par(18, Auswahl);
    d4    = Par(19, Auswahl);
    nu0   = Par(20, Auswahl);

    MuMin = 0;
    MuMax = 10;
    mu_paramVec = MuMin:SW:MuMax;
    nMu   = length(mu_paramVec);

    ebetaStr   = strrep(sprintf('_ebeta%2.3f', ebeta), '.', 'dot');
    muStr      = sprintf('_mu%d_mu%d', round(MuMin), round(MuMax));
    auswahlStr = sprintf('_Auswahl%d', Auswahl);
    filename   = ['Workspace', ebetaStr, muStr, auswahlStr, '.mat'];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Pruefen ob Berechnung noetig

    needsCompute = true;
    if exist(filename, 'file') == 2 && loadMat
        fileVars = who('-file', filename);
        if all(ismember({'CharExRe1Cor','CharExRe2Cor','CharExIm1Cor','CharExIm2Cor', ...
                         'CharExRe','CharExIm','tableCharEx'}, fileVars))
            load(filename);
            needsCompute = false;
            fprintf('  Loaded from file.\n');
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Berechnung

    if needsCompute

        idx        = 1;
        n          = 1;
        cntN       = 1;
        prevAbsDiffCharMult = 0;

        Diagonal    = diag(ones(AnzGl,1));
        Monodromie  = zeros(AnzGl);
        CharMult    = zeros(nMu, AnzGl+1);
        CharExRe    = zeros(nMu, AnzGl);
        CharExIm    = zeros(nMu, AnzGl);
        CharExImRaw = zeros(nMu, AnzGl);
        CharEx      = zeros(nMu, AnzGl*6);
        cond_A      = nan(nMu, 1);
        buffer.Pos  = zeros(nPairs, 1);

        [~,A] = SchlagDGL(0, 0, gamma, d2, d3, d4, 0, ebeta, nu0, Blatt);

        [freq, ~] = sort(imag(eig(A)).', 'descend'); %#ok<UDIM>
        damp = real(eig(A));

        if freq(1) - fix(freq(1)) < 0.5
            freqInt = round(freq);
        else
            freqInt = -round(freq);
        end

        if Blatt == 5
            freqInt = freqInt([1,2,3,6,7,4,5,8,9,10]);
        end

        options2 = odeset('RelTol',1e-10,'AbsTol',1e-12,'MaxStep',1e-3);

        for mu_param = mu_paramVec

            if konstant == 1
                for k = 1:AnzGl
                    sol = ode45(@(psi,x) SchlagDGLkonstant(psi,x,gamma,d2,d3,d4,mu_param,ebeta,nu0,Blatt), ...
                        [t0,T], Diagonal(:,k), options2);
                    Monodromie(:,k) = deval(sol,T);
                end
            else
                for k = 1:AnzGl
                    sol = ode45(@(psi,x) SchlagDGL(psi,x,gamma,d2,d3,d4,mu_param,ebeta,nu0,Blatt), ...
                        [t0,T], Diagonal(:,k), options2);
                    Monodromie(:,k) = deval(sol,T);
                end
                cond_A(idx) = cond(Monodromie);
            end

            % Charakteristische Multiplikatoren
            charMult     = (eig(Monodromie))';
            [~,idxSort]  = sort(real(charMult));
            charMultSort = charMult(idxSort);

            CharMult(idx,:) = [charMultSort, mu_param];

            if mu_param == 0
                Im0 = 1/T * angle(charMultSort);
            end

            % Charakteristische Exponenten
            Eig.Real = 1/T * log(abs(charMultSort));
            Eig.Imag = 1/T * atan(imag(charMultSort) ./ real(charMultSort));

            % correctImagValues paarweise
            ImagCorrected    = zeros(1, nPairs);
            ImagCorrectedNeg = zeros(1, nPairs);

            for idxC = 1:nPairs
                idxVec = 2*idxC-1 : 2*idxC;
                EigTemp.Real   = Eig.Real(idxVec);
                EigTemp.Imag   = Eig.Imag(idxVec);
                bufferTemp.Pos = buffer.Pos(idxC);

                [EigTemp, bufferTemp] = correctImagValues(EigTemp, bufferTemp);

                Eig.Real(idxVec)     = EigTemp.Real;
                Eig.Imag(idxVec)     = EigTemp.Imag;
                Eig.ImagSort(idxVec) = EigTemp.ImagSort;
                ImagCorrected(idxC)    = EigTemp.ImagCorrected;
                ImagCorrectedNeg(idxC) = EigTemp.ImagCorrectedNeg;
                buffer.Pos(idxC) = bufferTemp.Pos;
            end

            Eig.ImagCorrected    = ImagCorrected;
            Eig.ImagCorrectedNeg = ImagCorrectedNeg;

            % Additionsterm — nAddVec muss Groesse 1 x AnzGl haben
            absDiffCharMult = abs(CharMult(idx,1)) > abs(CharMult(idx,b1));
            if idx > 1 && abs(prevAbsDiffCharMult - absDiffCharMult) > 0.1
                n    = n + 0.5;
                cntN = cntN + 1;
            end
            prevAbsDiffCharMult = absDiffCharMult;

            nAdd     = n * 2*pi / T;
            nAddVec  = [-nAdd*ones(1,nPairs), nAdd*ones(1,nPairs)];  % 1 x AnzGl

            ImagEigSortN = [Eig.ImagCorrectedNeg, Eig.ImagCorrected] + nAddVec;

            CharEx(idx,:)     = [Eig.Real, Eig.Imag, Eig.ImagSort, ImagEigSortN, ...
                                  real(charMultSort), imag(charMultSort)];
            CharExRe(idx,:)    = Eig.Real;
            CharExIm(idx,:)    = Eig.Imag + nAddVec;
            CharExImRaw(idx,:) = Eig.Imag;

            idx = idx + 1;
        end

        %% Tabelle
        lEig = length(Eig.Real);
        cellRe   = regexp(sprintf('RealCharExp%02d,',  1:lEig), ',', 'split');
        cellIm   = regexp(sprintf('ImagCharExp%02d,',  1:lEig), ',', 'split');
        cellImS  = regexp(sprintf('ImagSort%02d,',     1:lEig), ',', 'split');
        cellImSN = regexp(sprintf('ImagSortN%02d,',    1:lEig), ',', 'split');
        cellReM  = regexp(sprintf('RealCharMult%02d,', 1:lEig), ',', 'split');
        cellImM  = regexp(sprintf('ImagCharMult%02d,', 1:lEig), ',', 'split');

        VarNames = ['muChar', cellRe(1:lEig), cellIm(1:lEig), cellImS(1:lEig), ...
                    cellImSN(1:lEig), cellReM(1:lEig), cellImM(1:lEig)];

        tableCharEx    = array2table([CharMult(:,end), CharEx], 'VariableNames', VarNames);
        PrintNames     = VarNames(contains(VarNames,'Char'));
        tableCharPrint = tableCharEx(:, PrintNames);

        %% Korrektur Realteile
        CharExRe1 = CharExRe(:,1);
        CharExRe2 = CharExRe(:,Blatt+1);

        CharExRePos      = max(CharExRe1, CharExRe2);
        CharExRe1_NegIdx = abs(CharExRePos - CharExRe1) > eps;
        offset           = CharExRe1(1);
        CharExReNeg      = -(CharExRePos - offset) + offset;

        CharExRe1Cor = CharExRe1;
        CharExRe1Cor(CharExRe1_NegIdx)  = CharExReNeg(CharExRe1_NegIdx);
        CharExRe2Cor = CharExRe2;
        CharExRe2Cor(~CharExRe1_NegIdx) = CharExReNeg(~CharExRe1_NegIdx);

        tableCharPrint.RealCharExp1Corrected = CharExRe1Cor;
        tableCharPrint.RealCharExp2Corrected = CharExRe2Cor;

        %% Korrektur Imaginaerteile
        CharExIm1 = CharExIm(:,1);
        CharExIm2 = CharExIm(:,Blatt+1);

        CharExIm1Cor = CharExIm1;
        CharExIm2Cor = CharExIm2;

        for k = 1:nMu
            if CharExRe1Cor(k) == CharExRe1(k)
                CharExIm1Cor(k) =  CharExIm1(k);
                CharExIm2Cor(k) = -CharExIm1(k);
            else
                CharExIm1Cor(k) = -CharExIm2(k);
                CharExIm2Cor(k) =  CharExIm2(k);
            end
        end

        tableCharPrint.ImagCharExp1Corrected = CharExIm1Cor;
        tableCharPrint.ImagCharExp2Corrected = CharExIm2Cor;

        %% Speichern
        excelfilename  = strrep(filename, '.mat', '.xlsx');
        excelfilename  = strrep(excelfilename, 'Workspace', 'CharactExponentenMultiplikatoren');
        excelfilename1 = fullfile(excelDir, excelfilename);

        save(filename, 'damp','freq','MuMin','SW','MuMax','mu_paramVec', ...
            'CharExRe','CharExIm','CharExImRaw','CharEx','CharMult', ...
            'nu0','Blatt','tableCharEx','tableCharPrint', ...
            'CharExRe1Cor','CharExRe2Cor','CharExIm1Cor','CharExIm2Cor');
        save('Workspace.mat', 'damp','freq','MuMin','SW','MuMax','mu_paramVec', ...
            'CharExRe','CharExIm','CharExImRaw','CharEx','CharMult', ...
            'nu0','Blatt','tableCharEx','tableCharPrint', ...
            'CharExRe1Cor','CharExRe2Cor','CharExIm1Cor','CharExIm2Cor');

        writetable(tableCharPrint, excelfilename1);

    end % needsCompute

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Abbildung

    rotorDescription = AuswahlInfo{Auswahl, 3};
    mu_vec = MuMin:SW:MuMax;
    cl = lines;

    figure(Auswahl); clf;
    pos0 = get(0,"defaultFigurePosition");
    set(gcf, "Position", [pos0(1:2), 1.1*pos0(3),pos0(4)] )
    tl = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
    title(tl, sprintf('Auswahl: %d; Blatt: %d; %s', Auswahl, Blatt, rotorDescription));

    ax2(1) = nexttile;
    plot(mu_vec, CharExRe1Cor,   '-',   'Color', cl(1,:), 'LineWidth', 1.5); hold on;
    plot(mu_vec, CharExRe2Cor,   '-',   'Color', cl(2,:), 'LineWidth', 1.5);
    plot(mu_vec, CharExRe(:,1),  'k--', 'LineWidth', 1);
    plot(mu_vec, CharExRe(:,b1), 'k-.', 'LineWidth', 1);
    grid on; ylabel('Real part (-)');
    legend('Re(Exp1) cor', sprintf('Re(Exp%d) cor', Blatt), ...
           'Re(Exp1) raw',  sprintf('Re(Exp%d) raw', Blatt), ...
           'Location','northeastoutside');

    ax2(2) = nexttile;
    plot(mu_vec, CharExIm1Cor,   '-',   'Color', cl(1,:), 'LineWidth', 1.5); hold on;
    plot(mu_vec, CharExIm2Cor,   '-',   'Color', cl(2,:), 'LineWidth', 1.5);
    plot(mu_vec, CharExIm(:,1),  'k--', 'LineWidth', 1);
    plot(mu_vec, CharExIm(:,b1), 'k-.', 'LineWidth', 1);
    grid on;
    xlabel('Advance ratio \mu (-)'); ylabel('Imaginary part (-)');
    legend('Im(Exp1) cor', sprintf('Im(Exp%d) cor', Blatt), ...
           'Im(Exp1) raw',  sprintf('Im(Exp%d) raw', Blatt), ...
           'Location','northeastoutside');

    linkaxes(ax2,'x');

end % Auswahl loop

if plotAll == 0
    return;
end

%% Diagnostik-Plots
figure(100);
ax1(1) = subplot(2,1,1);
plot(MuMin:SW:MuMax, CharExRe(:,1),     '*-', ...
     MuMin:SW:MuMax, CharExRe(:,b1),    'o-', ...
     MuMin:SW:MuMax, CharExIm(:,1),     '*-', ...
     MuMin:SW:MuMax, CharExIm(:,b1),    'o-', ...
     MuMin:SW:MuMax, CharExImRaw(:,1),  '*-', ...
     MuMin:SW:MuMax, CharExImRaw(:,b1), 'o-');
legend('Re(Exp1)',sprintf('Re(Exp%d)',b1),'Im(Exp1)',sprintf('Im(Exp%d)',b1), ...
       'ImRaw(Exp1)',sprintf('ImRaw(Exp%d)',b1),'Location','SouthWest');
grid on;
%title(sprintf('Auswahl: %d; Blatt: %d', Auswahl, Blatt));
title(tl, {sprintf('Auswahl: %d; Blatt: %d', Auswahl, Blatt), rotorDescription},'FontSize', 9);

ax1(2) = subplot(2,1,2);
plot(MuMin:SW:MuMax, CharExRe1Cor, '*-', ...
     MuMin:SW:MuMax, CharExRe2Cor, 'o-', ...
     MuMin:SW:MuMax, CharExIm1Cor, '*-', ...
     MuMin:SW:MuMax, CharExIm2Cor, 'o-');
grid on;
linkaxes(ax1,'x');

figure(102);
lineCell = {'-','--','-.','-','--','-.'};
for idx = 1:Blatt
    if Blatt > 1; subplot(Blatt,1,idx); end
    plot(MuMin:SW:MuMax, CharExRe(:,idx),      lineCell{idx}, ...
         MuMin:SW:MuMax, CharExRe(:,idx+Blatt), lineCell{idx+Blatt});
    grid on;
    legend(cellRe([idx, idx+Blatt]));
    if idx == 1
        title(sprintf('Auswahl: %d; Blatt: %d', Auswahl, Blatt));
    end
end
xlabel('\mu (-)');