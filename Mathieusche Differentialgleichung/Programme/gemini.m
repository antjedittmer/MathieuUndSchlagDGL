%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation of the Strutt Diagrams within the bounds of nu_02 and nu_C2
% Combined script including all required functions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;

% --- Parameters and Setup ---
loadMat = 0;  % load mat-file if results with same D already exist
SW = 0.1; % step width
unt0 = 0; untC = 0;
ob0 = 9; obC = 9;
fDir = 'figureFolder1'; % Folder for figures
dDir = 'dataFolder'; % Folder for mat-files
excelDir = 'dataFolder1';

% Create necessary directories
if ~isdir(fDir), mkdir(fDir); end
if ~isdir(dDir), mkdir(dDir); end
if ~isdir(excelDir), mkdir(excelDir); end

% Number of equations in the ODE system
Nz = 2;
% Parameters
DVec = 0.15; % Damping parameter D
t0 = 0.0;
T = 2*pi;

for dIdx = 1: length(DVec)
    D = DVec(dIdx);
    
    % Initial conditions for Monodromy matrix calculation
    Diagonal = diag(ones(Nz,1));
    
    % Matlab mat-file name
    matName = [strrep(sprintf('STRUTTscheKarte_D%2.1e_SW%2.1e',D,SW),'.','dot'),...
        '.mat'];
    if unt0 == 0
        matName = strrep(matName,'.mat','_unt0.mat');
    end
    fileName = fullfile(dDir,matName);
    
    % --- Preparation of arrays ---
    nu_02_values = unt0:SW:ob0;
    nu_C2_values = untC:SW:obC;
    lenNu02 = length(nu_02_values); 
    lenNuC2 = length(nu_C2_values); 
    lenNu = lenNu02 * lenNuC2;
    
    lenNuDiag = min(lenNu02,lenNuC2); % Number of nu_02 == nu_C2 values
    CharEx = zeros(lenNuDiag,Nz*4+2); % Stores [nu02, nuC2, Re1, Re2, Im_min, Im_max, Im_corr_N1, Im_corr_N2, mu1, mu2]
    nAddVector = nan(lenNuDiag,1);
    plotwertstabil = zeros(lenNu,3); % Stores [nu02, nuC2, stability_flag]
    
    nuCSwitchVec = [0.2,1,1.5^2,2^2,2.5^2] - 0.15; % Switch points for n
    n = 0; % Initial n for correction
    cntN = 1;
    buffer.Pos = 0; % Buffer for continuous positive Imaginary part
    buffer.Neg = 0; % Buffer for continuous negative Imaginary part
    
    lidx = 1; % Index for stability points (plotwertstabil)
    oidx = 1; % Index for diagonal points (CharEx)
    
    
    % --- Main Calculation Loop ---
    if exist(fileName,'file') == 2 && loadMat == 1
        load(fileName,'Char*','plotwert*', 'nu*','nAddVector');
    else
        options = odeset('RelTol',1e-10,'AbsTol',1e-12);
        
        for nu_02 = nu_02_values
            for nu_C2  = nu_C2_values
                % Numerical integration to find Monodromy matrix
                Monodromie = zeros(Nz);
                for k = 1 : Nz
                    % Integrate: d/dpsi [phi; dphi/dpsi] = f([phi; dphi/dpsi], psi)
                    sol = ode45(@(psi,x)MathieuDGL(x,psi,D,nu_02,nu_C2),[t0,T],Diagonal(:,k),options);
                    MonoVek  = deval(sol,T);
                    Monodromie(:,k) = MonoVek;
                end
                
                % Characteristic multipliers (eigenvalues of M)
                eP = eig(Monodromie);
                
                %% Characteristic exponents (only calculated on the diagonal nu_02 = nu_C2)
                if abs(nu_C2 - nu_02) < SW/2 % Check for diagonal
                    
                    % Calculate real and imaginary parts
                    Eig.Real = 1/T * log(abs(eP));
                    Eig.Imag = 1/T * atan2(imag(eP),real(eP)); % Use atan2 for full circle
                    Eig.ImagSort = sort(Eig.Imag);
                    
                    % Correct imaginary part for a continuous curve
                    [ImagEigCorrected, ImagEigCorrectedNeg, buffer] = correctImagValues(Eig.ImagSort, buffer);
                    
                    % Apply n*2*pi/T correction
                    if cntN <= length(nuCSwitchVec) && ... 
                            nu_C2 > nuCSwitchVec(cntN) && ... 
                            abs(Eig.Real(1) - Eig.Real(2)) > eps % Real parts must be equal at boundary
                        n =  n + 0.5;
                        cntN = cntN+1;
                    end
                    nAdd = n*2*pi/T; 
                    ImagEigSortN = [ImagEigCorrectedNeg, ImagEigCorrected] + [-nAdd,nAdd];
                    
                    CharEx(oidx,:) = [nu_02, nu_C2, Eig.Real', ImagEigCorrectedNeg, ImagEigCorrected, ImagEigSortN, eP'];
                    nAddVector(oidx) = nAdd;
                    oidx = oidx + 1;
                end
                
                %% Stability check: max(|mu|) < 1
                if  max(abs(eP)) < 1
                    b = 1; % Stable
                    plotwertstabil(lidx,:) = [nu_02,nu_C2,b];
                end
                lidx = lidx + 1;
            end
        end
        % Delete purely 0-rows and save
        plotwertstabil =  plotwertstabil(any(plotwertstabil,2),:);
        CharEx = CharEx(any(CharEx,2),:);
        nAddVector = nAddVector(~isnan(nAddVector));
        save(fileName,'Char*','plotwert*', 'nu*','nAddVector');
    end
    
    % --- Data Output to Excel ---
    try
        CharExTable = array2table(CharEx(:,[1:4,7,8]),'VariableNames',...
            {'nu02','nu_C2','Eig.Real1','Eig.Real2', 'ImagEig1_corr', 'ImagEig2_corr'});
        CharExTable.RealCharExp1 = real(CharEx(:,9));
        CharExTable.ImagCharExp1 = imag(CharEx(:,9));
        CharExTable.RealCharExp2 = real(CharEx(:,10));
        CharExTable.ImagCharExp2 = imag(CharEx(:,10));
        
        excelfilename = strrep(fileName,'.mat','CharExAll.xlsx');
        excelfilename1 = fullfile(excelDir, strrep(excelfilename,dDir,''));
        writetable(CharExTable,excelfilename1)
    catch ME
        disp(['Error writing Excel file: ', ME.message]);
    end
    
    % --- Graphical representation ---
    cl = lines;
    fs = 10.5; %Fontsize
    
    % Figure 1: Strutt Diagram and Exponents
    hf = figure(dIdx);
    hf.Position = [10 10 600 600];
    tiledlayout(11,1);
    
    % 1. Strutt diagram
    h(1) = nexttile([7 1]);
    scatter(plotwertstabil(:,1),plotwertstabil(:,2),5,'filled')
    ylabel('$\rm{Parameter} \; \nu_C^2 \;\;\rm{[-]}$','interpreter','latex','Position', [-0.5 obC/2],'FontSize', fs+2);
    grid on;
    title(['Strutt Diagram for D = ', num2str(D)], 'Interpreter', 'latex', 'FontSize', fs+4);
    
    % 2. Real parts of characteristic exponents
    h(2) = nexttile([2 1]);
    xachse = CharEx(:,1); 
    plot(xachse,CharEx(:,3:4))
    ylabel('$Re(s_R) \;\; \rm{[-]}$','interpreter','latex','FontSize', fs); 
    grid on;
    
    % 3. Imaginary parts of characteristic exponents (corrected)
    h(3) = nexttile([2 1]);
    plot(xachse,CharEx(:,7:8)); hold on;
    grid on;
    ylim(h(3),[-0.25,3]);
    xlabel('$\rm{Parameter} \; \nu_0^2 \;\;\rm{[-]}$','interpreter','latex','FontSize', fs+2);
    ylabel('$Im(s_R) \;\; \rm{[-]}$','interpreter','latex','FontSize', fs); 
    
    linkaxes(h,'x')
    for idxH = 1: length(h)
        set(h(idxH),'TickLabelInterpreter','Latex','FontSize',fs)
    end
    pngname = fullfile(fDir,strrep(matName,'.mat','_Strutt'));
    print(pngname, '-dpng')
    
    % Figure 2: Characteristic Exponents Detail
    hf1 = figure(dIdx+100);
    tiledlayout(2,1); 
    
    % Real part
    h(1) = nexttile; 
    plot(xachse,CharEx(:,3:4)); grid on;
    title({['$\ddot{\phi} + 2D \dot{\phi} + (\nu^2_0 + \nu^2_C \cos(\psi))\phi = 0$, D = ', num2str(D,2), ', $\nu_C = \nu_0 $']...
        '$\Re(s_{Ri}) = 1/2\pi\ln(|\mu_{Ri}|)$'},'interpreter','latex','FontSize', fs+2)
    ylabel('$\Re(s_R) \;\; \rm{[-]}$','interpreter','latex','FontSize', fs+2);
    
    % Imaginary part with n-jumps
    h(2) = nexttile; 
    plot(xachse,CharEx(:,7:8)); hold on;  grid on;
    try
        % Find where nAddVector increases (where n is incremented)
        idxMVec = [0;diff(nAddVector)]>0;
        idxM = xachse(idxMVec);
        for idxPl = 1: length(idxM)
            plot(idxM(idxPl)*ones(2,1), get(h(2),'ylim'),'--','Color',0.5*ones(3,1),'LineWidth',1);
        end
        
        legend('$s_{R1}$','$s_{R2}$','increase in n', 'interpreter','latex','Location','SouthEast','FontSize', fs+1)
        title('$\Im(s_{Ri}) = 1/2\pi (\arctan(\Im(\mu_{Ri})/\Re(\mu_{Ri})) + n \cdot 2\pi$','interpreter','latex','FontSize', fs+2)
    catch
    end
    xlabel('$\rm{Parameter} \; \nu_0^2 = \nu_C^2 \;\;\rm{[-]}$','interpreter','latex','FontSize', fs+2);
    ylabel('$\Im(s_R) \;\; \rm{[-]}$','interpreter','latex','FontSize', fs); 
    pngname = fullfile(fDir,strrep(strrep(strrep(matName,'.mat',''),'STRUTTscheKarte','CharExp'),'_unt0',''));
    print(pngname, '-dpng')
end

%% --- Helper Functions ---

% Function to convert the second-order ODE to two first-order ODEs
function dxdpsi = MathieuDGL(x, psi, D, nu_02, nu_C2)
% x(1) = phi
% x(2) = dphi/dpsi
% dphi/dpsi = x(2)
% d2phi/dpsi2 = -2*D*x(2) - (nu_0^2 + nu_C^2*cos(psi))*x(1)

dxdpsi = zeros(2,1);
dxdpsi(1) = x(2);
dxdpsi(2) = -2*D*x(2) - (nu_02 + nu_C2*cos(psi))*x(1);
end

function  [ImagEigCorrected, ImagEigCorrectedNeg,buffer,Eig] = correctImagValuesEig(eP,buffer,T)
%[ImagEigSortN,ImagEigCorrected, ImagEigCorrectedNeg, buffer.Pos,buffer.Neg,cntN] = correctImagValues(tmp,tmpNeg,buffer,buffer.Neg,cntN,nuCSwitchVec,n,nu_C2,T,RealEig)
% correctImagValues korrigiert die Werte

Eig.Real = 1/T * log(abs(eP));
Eig.Imag = 1/T * atan(imag(eP)./real(eP));
Eig.ImagAngle = 1/T * angle(imag(eP)./real(eP));
Eig.ImagSort = sort(Eig.Imag);

% Imaginaeranteil kontinuierlich steigend oder fallend
tmp = Eig.ImagSort(2);
tmpNeg = Eig.ImagSort(1);


if buffer.Pos <= tmp || (abs(tmp) < 10^-5) % Wert uebernehmen
    ImagEigCorrected = tmp; % steigender pos. Wert
    ImagEigCorrectedNeg = tmpNeg;  % fallender neg. Wert
    buffer.Pos = 0;
    buffer.Neg = 0;
else  % 'korrigierten' Wert nehmen fuer kontinuierlichen Verlauf
    % Wird der 'else'-Zweig getriggert, ist der buffer
    % bei 90 Grad (1.57 rad). Der positive Imaginaerteil ist damit
    % die Summe aus 180 Grad und dem negativen Winkel, dessen
    % Wert von -90 Grad zu 0 Grad laeuft.
    ImagEigCorrected = 2*buffer.Pos + tmpNeg; % korrigierter pos. Wert
    ImagEigCorrectedNeg = 2*buffer.Neg + tmp; % korrigierter neg. Wert
    % vecSwitch = 1; %vecSwitch(oidx)= 1;

end

% buffer.Pos ueberschreiben mit aktuellem Wert
buffer.Pos = max(tmp,buffer.Pos); % Maximum speichern
buffer.Neg = min(tmpNeg,buffer.Neg); % Minimum speichern

% % Additionsterm n*2*pi/T
% if cntN <= length(nuCSwitchVec) && ... % Sicherheitscheck, damit n nicht groesser wird als die Anzahl der 'Switchstellen'
%         nu_C2 > nuCSwitchVec(cntN) && ... % n nur bei erreichen der Switchstellen umstellen
%         abs(RealEig(1) - RealEig(2))> eps % die Realteile muessen gleich sein
%     n =  n + 0.5;
%     cntN = cntN+1;
% end
%
% nAdd = n*2*pi/T;
% ImagEigSortN = [ImagEigCorrectedNeg, ImagEigCorrected] + [-nAdd,nAdd];
