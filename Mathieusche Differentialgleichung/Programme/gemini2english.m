%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation of Strutt's Charts within the bounds of nu_02 and nu_C2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;
loadMat = 1;  % Load mat-file if results with the same D exist
SW = 0.1; % stepwidth
unt0 = 0; % Lower bound for nu_02
untC = 0; % Lower bound for nu_C2
ob0 = 9;  % Upper bound for nu_02
obC = 9;  % Upper bound for nu_C2

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
DVec = 0.15; % Damping coefficient vector
t0 = 0.0; % Initial time
T = 2*pi; % Period of excitation

for dIdx = 1: length(DVec)
    D = DVec(dIdx); % Current damping coefficient
    
    % Initial conditions (Identity matrix columns for Monodromy)
    Diagonal = diag(ones(Nz,1));
    
    % Matlab mat-file Name
    matName = [strrep(sprintf('STRUTT_Chart_D%2.1e_SW%2.1e',D,SW),'.','dot'),...
        '.mat'];
    if unt0 == 0
        matName = strrep(matName,'.mat','_unt0.mat');
    end
    fileName = fullfile(dDir,matName);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Preparation of arrays to be filled in the loops for performance
    % acceleration
    lenNu02 = length(unt0:SW:ob0); % Number of Nu_02 values
    lenNuC2 = length(untC:SW:obC); % Number of Nu_C2 values
    lenNu = lenNu02 * lenNuC2; % Total number of Nu_02/Nu_C2 combinations
    lenNuDiag = min(lenNu02,lenNuC2); % Number of values where nu_02 == nu_C2 (Diagonal)
    
    Monodromie = zeros(Nz);
    % nu, nc, Real1, Real2, ImagSortMin, ImagSortMax, ImagCorrectedNeg, ImagCorrectedPos, eP1, eP2
    CharEx = zeros(lenNuDiag,Nz*4+2); 
    nAddVector = nan(lenNuDiag,1);
    plotwertstabil = zeros(lenNu,3); % [nu_02, nu_C2, stability_flag]
    
    % Switch points for the imaginary part correction (multiples of pi/T)
    nuCSwitchVec = [0.2,1,1.5^2,2^2,2.5^2] - 0.15;
    n = 0;
    cntN = 1;
    buffer.Pos = 0;
    buffer.Neg = 0; 
    vecSwitch = zeros(lenNu02,1);
    noF = 1; % Figure number
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if exist(fileName,'file') == 2 && loadMat == 1
        load(fileName,'Char*','plotwert*', 'nu*','nAddVector');
    else
        % Counter variables
        lidx = 1; % Index for characteristic multipliers (stability flag)
        oidx = 1; % Index for characteristic exponents (diagonal plot)
        
        for nu_02 = unt0:SW:ob0
            for nu_C2  = untC:SW:obC
                % The loop numerically integrates Nz times for the column
                % vectors of the identity matrix. The solution at T is 
                % written to MonoVek, which is then assembled into the
                % Monodromy matrix.
                
                options = odeset('RelTol',1e-10,'AbsTol',1e-12);
                for k = 1 : Nz
                    % The system of ODEs: MathieuDGL(psi,x,D,nu_02,nu_C2)
                    sol = ode45(@(psi,x)MathieuDGL(psi,x,D,nu_02,nu_C2),[t0,T],Diagonal(:,k),options);
                    MonoVek  = deval(sol,T);
                    Monodromie(:,k) = MonoVek;
                end
                
                % Characteristic Multipliers (Eigenvalues of the Monodromy matrix)
                eP = eig(Monodromie);
                
                %% Characteristic Exponents
                if nu_C2 == nu_02 % Only calculated along the diagonal
                    % Calculate Real and Imaginary parts of the exponents
                    Eig.Real = 1/T * log(abs(eP));
                    Eig.Imag = 1/T * atan(imag(eP)./real(eP));
                    
                    % Correct the Imaginary part for a continuous curve
                    [Eig,buffer] = correctImagValues(Eig,buffer);
                    
                    % % Addition term n*2*pi/T
                    if cntN <= length(nuCSwitchVec) && ... % Safety check: n shouldn't exceed the number of 'switch points'
                            nu_C2 > nuCSwitchVec(cntN) && ... % Only adjust n when reaching a switch point
                            abs(Eig.Real(1) - Eig.Real(2))> eps % Imaginary parts must be equal (at the switch point)
                        n =  n + 0.5;
                        cntN = cntN+1;
                    end
                    nAdd = n*2*pi/T; % The term added to ensure continuity
                    
                    % Apply the addition term
                    ImagEigSortN = [Eig.ImagCorrectedNeg, Eig.ImagCorrected] + [-nAdd,nAdd];
                    
                    % Store the results: [nu02, nuC2, Real1, Real2, ImagSortMin, ImagSortMax, ImagCorrectedNeg+Term, ImagCorrectedPos+Term, eP1, eP2]
                    CharEx(oidx,:) = [nu_02, nu_C2, Eig.Real', min(Eig.ImagSort), max(Eig.ImagSort), ImagEigSortN, eP'];
                    nAddVector(oidx) = nAdd;
                    oidx = oidx + 1;
                end
                
                %% Stability flag: 1 for stable combinations of nu_02 and nu_C2 based on characteristic multipliers
                if  max(abs(eP)) < 1 % Stable condition: all multipliers have magnitude < 1
                    b = 1;
                    plotwertstabil(lidx,:) = [nu_02,nu_C2,b];
                end
                lidx = lidx + 1;
            end
        end
        % Delete purely zero rows (from initialization)
        plotwertstabil =  plotwertstabil(any(plotwertstabil,2),:);
        save(fileName,'Char*','plotwert*', 'nu*','nAddVector');
    end
    
    %% Export to Excel
    try
        % Create table with characteristic exponent data
        CharExTable = array2table(CharEx(:,[1:4,7,8]),'VariableNames',...
            {'nu02','nu_C2','Eig.Real1','Eig.Real2', 'ImagEig1', 'ImagEig2'});
        
        % Append Real and Imaginary parts of the multipliers (eP)
        CharExTable.RealCharExp1 = real(CharEx(:,9));
        CharExTable.ImagCharExp1 = imag(CharEx(:,9));
        CharExTable.RealCharExp2 = real(CharEx(:,10));
        CharExTable.ImagCharExp2 = imag(CharEx(:,10));
        
        excelfilename = strrep(fileName,'.mat','CharExAll.xlsx');
        excelfilename1 = strrep(excelfilename,dDir,excelDir);
        writetable(CharExTable,excelfilename1)
    catch
        warning('Failed to write characteristic exponent data to Excel.');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Graphical Representation: Strutt's Chart (Plotting stable nu_02/nu_C2 pairs) and characteristic exponent plots
    cl = lines;
    fs = 10.5; % Fontsize
    hf = figure(dIdx);
    hf.Position = [10 10 600 600];
    
    % Use tiledlayout for modern MATLAB, fallback to subplot
    try tiledlayout(11,1); catch, end 
    
    % Strutt's Chart Plot (Top 7/11 parts)
    try h(1) = nexttile([7 1]); catch, h(1) = subplot(11,1,[1,7]);  end 
    scatter(plotwertstabil(:,1),plotwertstabil(:,2),5,'filled')
    ylabel('$\rm{Parameter} \; \nu_C^2 \;\;\rm{[-]}$','interpreter','latex','Position', [-0.5 obC/2],'FontSize', fs+2);
    grid on;
    
    % Real Parts of Characteristic Exponents (Middle 2/11 parts)
    try h(2) = nexttile([2 1]); catch, h(2) = subplot(11,1,[8,9]); end
    xachse = CharEx(:,1); 
    plot(xachse,CharEx(:,3:4))
    ylabel('$Re(s_R) \;\; \rm{[-]}$','interpreter','latex','FontSize', fs); 
    grid on;
    
    % Imaginary Parts of Characteristic Exponents (Bottom 2/11 parts)
    try h(3) = nexttile([2 1]); catch, h(3) = subplot(11,1,[10,11]); end 
    plot(xachse,CharEx(:,7:8)); hold on;
    grid on;
    ylim(gca,[-0.25,3]);
    grid on;
    xlabel('$\rm{Parameter} \; \nu_0^2 \;\;\rm{[-]}$','interpreter','latex','FontSize', fs+2);
    ylabel('$Im(s_R) \;\; \rm{[-]}$','interpreter','latex','FontSize', fs); 
    
    % Link x-axes of the three plots
    linkaxes(h,'x')
    for idxH = 1: length(h)
        set(h(idxH),'TickLabelInterpreter','Latex','FontSize',fs)
    end
    
    % Save combined figure
    pngname = fullfile(fDir,strrep(matName,'.mat',''));
    print(pngname, '-dpng')
    
    % Separate figure for characteristic exponents
    hf1 = figure(dIdx+100);
    try tiledlayout(2,1); catch, end
    
    % Real part
    try h(1) = nexttile; catch, h(1) = subplot(2,1,1);  end 
    plot(xachse,CharEx(:,3:4)); grid on;
    title({['$\ddot{\phi} + 2D \dot{\phi} + (\nu^2_0 + \nu^2_C \cos(\psi))\phi = 0$, D = ', num2str(D,2), ', $\nu_C = \nu_0 $']...
        '$\Re(s_{Ri}) = 1/2\pi\ln(|\mu_{Ri}|)$, $\mu_{Ri}$: Eigenvalue monodromy matrix'},'interpreter','latex','FontSize', fs+2) 
    ylabel('$\Re(s_R) \;\; \rm{[-]}$','interpreter','latex','FontSize', fs+2);
    
    % Imaginary part
    try h(2) = nexttile; catch, h(2) = subplot(2,1,2);  end 
    plot(xachse,CharEx(:,7:8)); hold on;  grid on;
    
    try
        % Add markers for jump points (change in nAddVector)
        idxMVec = [0,diff(nAddVector)]>0;
        idxM = xachse(idxMVec);
        for idxPl = 1: length(idxM)
            plot(idxM(idxPl)*ones(2,1), get(gca,'ylim'),'Color',0.5*ones(3,1),'LineWidth',1);
        end
        for idxPl = 1: length(idxM)
            tmp = nAddVector(idxMVec);
            text( idxM(idxPl) + 0.1, max(get(gca,'ylim'))-1, num2str(tmp(idxPl),2),'interpreter','latex','FontSize', fs+1);
        end
        legend('$s_{R1}$','$s_{R2}$','increase in m', 'interpreter','latex','Location','SouthEast','FontSize', fs+1)
        title('$\Im(s_{Ri}) = 1/2\pi (\arctan(\Im(\mu_{Ri})/\Re(\mu_{Ri})) + m$','interpreter','latex','FontSize', fs+2)
    catch
         % Keep this catch block for serious unexpected errors, 
        warning('Failed to add markers and legend to the imaginary exponent plot.'); 
    end
    xlabel('$\rm{Parameter} \; \nu_C^2 \;\;\rm{[-]}$','interpreter','latex','FontSize', fs+2);
    ylabel('$\Im(s_R) \;\; \rm{[-]}$','interpreter','latex','FontSize', fs+2); 
    
    % Save characteristic exponent figure
    pngname = fullfile(fDir,strrep(strrep(strrep(matName,'.mat',''),'STRUTT_Chart','CharExp'),'_unt0',''));
    print(pngname, '-dpng')
end

%% =========================================================================
%  LOCAL FUNCTIONS
%  =========================================================================

function  [Eig, buffer] = correctImagValues(Eig, buffer)
% correctImagValues corrects the imaginary parts for continuous progression.
% 
% Inputs:
% - Eig: Struct with vector Eig.Imag containing the two imaginary parts.
% - buffer: Struct containing the last values for maximum (Pos) and minimum (Neg).
% 
% Outputs:
% - Eig: Struct with appended, corrected imaginary parts.
% - buffer: Struct overwritten with new max/min values.

% Two checks for the correct format
if nargin~= 2
    error('Two inputs are expected: The current imaginary part of the eigenvalues and the buffer with the last values');
end
if max(size(Eig.Imag))~= 2 || min(size(Eig.Imag))~= 1
    error('The current imaginary parts of an eigenvalue pair is expected');
end

% Sort imaginary parts
Eig.ImagSort = sort(Eig.Imag);

% Imaginary part continuously increasing or decreasing
tmp = Eig.ImagSort(2); % The positive/upper value
tmpNeg = Eig.ImagSort(1); % The negative/lower value

% Check if the current positive value (tmp) is greater than the previous maximum (buffer.Pos)
if buffer.Pos <= tmp || (abs(tmp) < 10^-5) % Accept the value
    Eig.ImagCorrected = tmp; % Increasing positive value
    Eig.ImagCorrectedNeg = tmpNeg;  % Decreasing negative value
else  % 'Correct' the value to maintain continuity
    % If the 'else' branch is triggered, the previous buffer was around
    % 90 degrees (pi/2 rad). The corrected positive imaginary part is the
    % sum of 180 degrees (2*buffer.Pos, approx pi) and the current negative angle.
    Eig.ImagCorrected = 2*buffer.Pos  + tmpNeg; % Corrected positive value
    Eig.ImagCorrectedNeg  = 2*buffer.Neg + tmp; % Corrected negative value
end

% Overwrite the buffer with the current value
buffer.Pos = max(tmp,buffer.Pos); % Save the Maximum
buffer.Neg = min(tmpNeg,buffer.Neg); % Save the Minimum
end

% -------------------------------------------------------------------------

function dxdt = MathieuDGL(psi, x, D, nu_02, nu_C2)
% MathieuDGL defines the system of ODEs for the damped Mathieu equation:
% phi'' + 2*D*phi' + (nu_0^2 + nu_C^2 * cos(psi))*phi = 0
%
% State vector x:
% x(1) = phi
% x(2) = dphi/dpsi (phi')

dxdt = zeros(2,1);

% x(1)' = x(2)
dxdt(1) = x(2);

% x(2)' = -2*D*x(2) - (nu_0^2 + nu_C^2 * cos(psi))*x(1)
dxdt(2) = -2*D*x(2) - (nu_02 + nu_C2 * cos(psi))*x(1);
end