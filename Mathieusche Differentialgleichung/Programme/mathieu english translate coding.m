%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation of the Strutt Diagrams within the bounds of nu_02 and nu_C2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;
loadMat = 1;  % load mat-file if results with same D already exist
SW = 0.5; % step width
unt0 = 0;
untC = 0;
ob0 = 9;
obC = 9;
fDir = 'figureFolder1'; % Folder for figures
if ~isdir(fDir) %#ok<ISDIR>
    mkdir(fDir)
end
dDir = 'dataFolder'; % Folder for mat-files
if ~isdir(dDir) %#ok<ISDIR>
    mkdir(dDir)
end
excelDir = 'dataFolder1';
if ~isdir(excelDir) %#ok<ISDIR>
    mkdir(excelDir);
end
% Number of equations in the ODE system
Nz = 2;
% Parameters
DVec = 0.15; %0.001; %[0.15; 0.001; 0.2]; 0.3; %
% strDVec = {sprintf('%2.2f',DVec(1)); sprintf('%2.3f',DVec(2)); ...
%     sprintf('%2.1f',DVec(3))};
t0 = 0.0;
T = 2*pi;
for dIdx = 1: length(DVec)
    D = DVec(dIdx);
    % Initial conditions
    Diagonal = diag(ones(Nz,1));
    % Matlab mat-file name
    matName = [strrep(sprintf('STRUTTscheKarte_D%2.1e_SW%2.1e',D,SW),'.','dot'),...
        '.mat'];
    if unt0 == 0
        matName = strrep(matName,'.mat','_unt0.mat');
    end
    fileName = fullfile(dDir,matName);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Preparation of arrays to be filled in loops for program acceleration
    lenNu02 = length(unt0:SW:ob0); % Number of Nu_02 values
    lenNuC2 = length(untC:SW:obC); % Number of Nu_C2 values
    lenNu = lenNu02 * lenNuC2; % Number of Nu_C2 combinations
    lenNuDiag = min(lenNu02,lenNuC2); % Number of nu_02 == nu_C2 values
    Monodromie = zeros(Nz);
    CharEx = zeros(lenNuDiag,Nz*4+2); % nu, nc, Real1,2, Imag12, Imag1,2, Pole
    nAddVector = nan(lenNuDiag,1);
    plotwertstabil = zeros(lenNu,3);
    %nuCSwitchVec = [0.1,1,1.5^2,2^2,2.5^2] - 0.1;
    nuCSwitchVec = [0.2,1,1.5^2,2^2,2.5^2] - 0.15;
    n = 0;
    cntN = 1;
    buffer.Pos = 0;
    bufferNeg = 0;
    vecSwitch = zeros(lenNu02,1);
    noF = 1; % Figure Number
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if exist(fileName,'file') == 2 && loadMat == 1
        load(fileName,'Char*','plotwert*', 'nu*','nAddVector');
    else
        % Counting variables
        lidx = 1; % Index for characteristic multipliers
        oidx = 1; % Index for characteristic exponents
        for nu_02 = unt0:SW:ob0
            for nu_C2  = untC:SW:obC
                % The loop numerically integrates Nz times one after another for the
                % column vectors of the identity matrix. The solution at T is
                % written into MonoVek and assembled into the Monodromy matrix
                options = odeset('RelTol',1e-10,'AbsTol',1e-12);
                for k = 1 : Nz
                    sol = ode45(@(psi,x)MathieuDGL(psi,x,D,nu_02,nu_C2),[t0,T],Diagonal(:,k),options);
                    MonoVek  = deval(sol,T);
                    Monodromie(:,k) = MonoVek;
                end
                % Characteristic multipliers (eigenvalues of the monodromy matrix)
                eP = eig(Monodromie);
                %% Characteristic exponents
                if nu_C2 == nu_02
                    % Calculate real and imaginary parts of the exponents
                    Eig.Real = 1/T * log(abs(eP));
                    Eig.Imag = 1/T * atan(imag(eP)./real(eP));
                    % Correct imaginary part for a continuous curve
                    [Eig,buffer] = correctImagValues(Eig,buffer);
                    % % Addition term n*2*pi/T
                    if cntN <= length(nuCSwitchVec) && ... % Safety check, so n does not become larger than the number of 'switch points'
                            nu_C2 > nuCSwitchVec(cntN) && ... % only change n when a switch point is reached
                            abs(Eig.Real(1) - Eig.Real(2)) > eps % the imaginary parts must be equal
                        n =  n + 0.5;
                        cntN = cntN+1;
                    end
                    nAdd = n*2*pi/T; %
                    ImagEigSortN = [Eig.ImagCorrectedNeg, Eig.ImagCorrected] + [-nAdd,nAdd];
                    CharEx(oidx,:) = [nu_02, nu_C2, Eig.Real', min(Eig.ImagSort), max(Eig.ImagSort), ImagEigSortN, eP'];
                    nAddVector(oidx) = nAdd;
                    oidx = oidx + 1;
                end
                %% 1 = stable combinations of nu_02 and nu_C2 based on the
                % characteristic multipliers
                if  max(abs(eP)) < 1 %eBetr(1)<1 && eBetr(2)<1
                    b = 1;
                    plotwertstabil(lidx,:) = [nu_02,nu_C2,b];
                end
                lidx = lidx + 1;
            end
        end
        % Delete purely 0-rows
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
    %% Graphical representation: Strutt diagram (Plotting all value pairs
    % nu_02/nu_C2 that lead to a stable solution) as well as the real parts of the
    % characteristic exponents
    cl = lines;
    fs = 10.5; %Fontsize
    hf = figure(dIdx);
    hf.Position = [10 10 600 600];
    try tiledlayout(11,1); catch, end% tiledlayout not available before R2007
    try h(1) = nexttile([7 1]); catch, h(1) = subplot(11,1,[1,7]);  end %
    % Strutt diagram
    scatter(plotwertstabil(:,1),plotwertstabil(:,2),5,'filled')
    ylabel('$\rm{Parameter} \; \nu_C^2 \;\;\rm{[-]}$','interpreter','latex','Position', [-0.5 obC/2],'FontSize', fs+2);
    grid on;
    % Real parts of characteristic exponents
    try h(2) = nexttile([2 1]); catch, h(2) = subplot(11,1,[8,9]); end
    xachse = CharEx(:,1); %unt0:SW:ob0;
    plot(xachse,CharEx(:,3:4))
    ylabel('$Re(s_R) \;\; \rm{[-]}$','interpreter','latex','FontSize', fs); %'Position', [-0.5 -D]
    grid on;
    try h(3) = nexttile([2 1]); catch, h(3) = subplot(11,1,[10,11]); end %nexttile([2 1]);
    plot(xachse,CharEx(:,7:8)); hold on;
    % plot(xachse,CharEx(:,5),'--', 'color',cl(1,:)) % For debugging
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
    try tiledlayout(2,1); catch, end% tiledlayout not available before R2007
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
    ylabel('$\Im(s_R) \;\; \rm{[-]}$','interpreter','latex','FontSize', fs); %'Position', [-0.5 -D]
    pngname = fullfile(fDir,strrep(strrep(strrep(matName,'.mat',''),'STRUTTscheKarte','CharExp'),'_unt0',''));
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