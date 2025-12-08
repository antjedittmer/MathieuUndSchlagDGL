%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stability Analysis (Strutt Chart) using Arnold's classic Floquet calculation
% logic (atan + manual phase correction) within the structure of the 
% Peters' code, plotting the Symmetric Frequency Deviation Curve.
% Differential Equation Form: phi'' + 2D*phi' + (nu_0^2 + nu_C^2*cos(psi))*phi = 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;
% --- Setup and Parameters ---
loadMat = 1; % Load mat-file if results with the same D are already available
SW = 0.1;    % Step width for nu_0^2 and nu_C^2 sweep
unt0 = 0;    % Lower bound for nu_0^2 (x-axis)
untC = 0;    % Lower bound for nu_C^2 (y-axis)
ob0 = 9;     % Upper bound for nu_0^2
obC = 9;     % Upper bound for nu_C^2
% Folder setup
fDir = 'figureFolder_Arnold_Classic_Symmetric';
if ~isfolder(fDir)
    mkdir(fDir)
end
dDir = 'dataFolder_Arnold_Classic_Symmetric';
if ~isfolder(dDir)
    mkdir(dDir)
end
excelDir = 'dataFolder_Arnold_Excel_Classic_Symmetric';
if ~isfolder(excelDir)
    mkdir(excelDir);
end
Nz = 2; % Number of equations in the DGL system (2 for a 2nd-order ODE)
DVec = 0.15; % Damping coefficient D
t0 = 0.0;
Omega = 1; % Normalized excitation frequency (Omega = 1)
T = 2*pi / Omega; % Period of the excitation (T = 2*pi)

for dIdx = 1: length(DVec)
    D = DVec(dIdx);
    
    % Initial conditions (Identity matrix for the State Transition Matrix)
    Diagonal = diag(ones(Nz,1));
    
    % MATLAB mat-file Name for all data
    matName = [strrep(sprintf('STRUTTscheKarte_Arnold_D%2.1e_SW%2.1e',D,SW),'.','dot'),...
        '.mat'];
    fileName = fullfile(dDir,matName);
    
    % --- Array Pre-allocation ---
    nu02_vals = unt0:SW:ob0;
    nuC2_vals = untC:SW:obC;
    lenNu02 = length(nu02_vals);
    lenNuC2 = length(nuC2_vals);
    lenNu = lenNu02 * lenNuC2;
    lenNuDiag = min(lenNu02,lenNuC2); % Number of diagonal points (nu_0^2 == nu_C^2)
    
    % Storage: nu02, nuC2, Re_s1, Re_s2, Im_s1/Omega, Im_s2/Omega
    CharEx = zeros(lenNuDiag, 6); 
    plotwertstabil = zeros(lenNu,3);
    
    lidx = 1; % Index for stability map (Strutt Chart)
    oidx = 1; % Index for characteristic exponents (diagonal only)

    % --- Arnold's Manual Tracking Initialization ---
    n = 0; % The correction term 'm' or 'n', initialized to zero to start from zero.
    buffer.Pos = 0;
    buffer.Neg = 0;
    
    if exist(fileName,'file') == 2 && loadMat == 1
        load(fileName,'CharEx','plotwertstabil');
    else
        disp(['Calculating Strutt Chart (Arnold/Symmetric Style) for D = ', num2str(D), '...']);
        for nu_02 = nu02_vals
            for nu_C2  = nuC2_vals
                
                options = odeset('RelTol',1e-10,'AbsTol',1e-12);
                
                % Solve the DGL system for the Monodromy Matrix columns
                for k = 1 : Nz        
                        sol = ode45(@(psi,x)MathieuDGL(psi,x,D,nu_02,nu_C2), [t0,T], Diagonal(:,k), options);
                        MonoVek  = deval(sol,T);
                        Monodromie(:,k) = MonoVek;
                end
                
                % Characteristic Multipliers (Eigenvalues of the Monodromie Matrix)
                eP = eig(Monodromie);
                
                % --- Characteristic Exponents (Only Diagonal nu_0^2 = nu_C^2) ---
                if abs(nu_C2 - nu_02) < SW/2 
                    
                    % 1. Arnold Style: Use atan and the Real part log(|mu|)
                    Eig.Real = 1/T * log(abs(eP)); % sigma = Re(s)
                    
                    % *** REPLACED atan2 WITH atan AS REQUESTED ***
                    Eig.Imag = 1/T * atan(imag(eP)./real(eP)); % raw Im(s) using atan ONLY
                    
                    % 2. Arnold Style: Correct Imaginary values for a continuous plot
                    [Eig,buffer] = correctImagValues(Eig,buffer);
                    
                    % 3. Arnold Style: Manual phase shift tracking (n-tracking)
                    % If the real parts become equal (resonance), increment n by 0.5.
                    if abs(Eig.Real(1) - Eig.Real(2)) < 1e-6 && nu_02 > 0.05 
                        n = n + 0.5;
                    end
                    nAdd = n*2*pi/T; % The phase shift term (m*2*pi/T)
                    
                    % Combine smoothed Im(s) and phase shift (nAdd)
                    % Eig.ImagCorrected is typically the positive frequency.
                    PhysFreq1 = Eig.ImagCorrected + nAdd; 
                    PhysFreq2 = Eig.ImagCorrectedNeg + nAdd; 

                    % Sort by Real part (stability measure)
                    [Eig_Re_Sort, idx_sort] = sort(Eig.Real); 
                    
                    % Store the least-damped frequency (s1) and the other (s2)
                    if Eig.Real(1) > Eig.Real(2)
                        % s1 is least damped
                        CharEx(oidx,:) = [nu_02, nu_C2, Eig.Real', PhysFreq1/Omega, PhysFreq2/Omega];
                    else
                        % s2 is least damped
                        CharEx(oidx,:) = [nu_02, nu_C2, Eig.Real', PhysFreq2/Omega, PhysFreq1/Omega];
                    end
                    oidx = oidx + 1;
                end
                
                % --- Stable Combinations (Strutt Chart) ---
                if  max(abs(eP)) < 1 % Stability condition: |mu| < 1
                    b = 1;
                    plotwertstabil(lidx,:) = [nu_02,nu_C2,b];
                end
                lidx = lidx + 1;
            end
        end
        % Delete pure zero rows
        plotwertstabil =  plotwertstabil(any(plotwertstabil,2),:);
        save(fileName,'CharEx','plotwertstabil');
    end
    
    % --- Optional: Excel Export ---
    try
        CharExTable = array2table(CharEx,'VariableNames',...
            {'nu02','nuC2','Re_s1','Re_s2', 'Freq_s1_norm', 'Freq_s2_norm'});
        excelfilename = strrep(fileName,'.mat','_CharExPhys.xlsx');
        excelfilename1 = strrep(excelfilename,dDir,excelDir);
        writetable(CharExTable,excelfilename1)
    catch ME
        warning(['Error during Excel export: ', ME.message]);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Graphical Representation (Symmetric Deviation Plot)
    cl = lines;
    fs = 10.5;
    hf = figure(dIdx);
    hf.Position = [10 10 600 700];
    
    % Data for plotting (Diagonal points)
    xachse = CharEx(:,1);
    
    h(1) = subplot(3,1,1);
    % Strutt Chart: nu_0^2 vs nu_C^2
    scatter(plotwertstabil(:,1),plotwertstabil(:,2),5,'filled');
    title({['Strutt Chart for $\ddot{\phi} + 2D \dot{\phi} + (\nu^2_0 + \nu^2_C \cos(\psi))\phi = 0$ ($D = $', num2str(D), ')']...
        '$\rm{Stability: } \max(|\mu|) < 1$'},'interpreter','latex','FontSize', fs+2);
    ylabel('Parameter $\nu_C^2$ $\rm{[-]}$','interpreter','latex','FontSize', fs+2);
    grid on;
    
    % Real parts of char. exponents: $\sigma = Re(s)$
    h(2) = subplot(3,1,2);
    plot(xachse,CharEx(:,3),'LineWidth', 1.5, 'Color', cl(1,:), 'DisplayName', '$\sigma_1$');
    hold on;
    plot(xachse,CharEx(:,4),'LineWidth', 1.5, 'Color', cl(2,:), 'DisplayName', '$\sigma_2$');
    yline(0, 'k--');
    title('Real Part (Damping) $\sigma = Re(s)$','interpreter','latex','FontSize', fs+2);
    ylabel('Damping Exponent $\sigma$ $\rm{[-]}$','interpreter','latex','FontSize', fs+2); 
    grid on;
    legend('Location','NorthEast', 'Interpreter','latex');

    % Imaginary parts: $\omega/\Omega = Im(s)/\Omega$
    h(3) = subplot(3,1,3); 
    
    % --- PLOTTING OF THE SYMMETRIC DEVIATION CURVE ---
    % CharEx(:,5) holds the normalized primary tracked frequency (omega1/Omega) from Arnold's logic.
    omega_primary_norm = abs(CharEx(:,5)); 
    
    % Calculate the absolute deviation from the principal resonance center (0.5)
    deviation_norm = abs(omega_primary_norm - 0.5);
    
    plot(xachse, deviation_norm,'LineWidth', 1.5, 'Color', cl(1,:), 'DisplayName', '$|\omega/\Omega - 0.5|$'); 
    hold on;
    grid on
    
    % Save the single frequency curve data to a new .mat file
    SingleFreqData.nu02_x_axis = xachse;
    SingleFreqData.deviation_normalized_frequency = deviation_norm;
    
    singleFreqFileName = fullfile(dDir,strrep(matName,'STRUTTscheKarte','DeviationFreq'));
    save(singleFreqFileName, 'SingleFreqData');
    
    disp(['Symmetric frequency deviation curve data saved to: ', singleFreqFileName]);
    % ---------------------------------------------------------------------
    
    % Reference line for the natural frequency deviation
    plot(xachse, abs(sqrt(xachse)/Omega - 0.5), 'k--', 'DisplayName', '$|\sqrt{\nu_0^2}/\Omega - 0.5|$');
    
    title('Normalized Frequency Deviation $|\omega/\Omega - 0.5|$ (Symmetric Curve)','interpreter','latex','FontSize', fs+2)
    xlabel('Parameter $\nu_0^2$ $\rm{[-]}$','interpreter','latex','FontSize', fs+2);
    ylabel('Deviation $\Delta\omega/\Omega$ $\rm{[-]}$','interpreter','latex','FontSize', fs+2); 
    legend('Location','NorthEast', 'Interpreter','latex');

    linkaxes(h,'x')
    for idxH = 1: length(h)
        set(h(idxH),'TickLabelInterpreter','Latex','FontSize',fs)
    end
    
    pngname = fullfile(fDir,strrep(matName,'.mat','_Exponents'));
    print(pngname, '-dpng')
end

% =========================================================================
% --- AUXILIARY FUNCTIONS (from the user's Arnold code snippet) ---
% The petersPhysicalFrequency function is removed.
% =========================================================================

% 1. Arnold's Manual Phase Correction Function
function  [Eig, buffer] = correctImagValues(Eig, buffer)
% correctImagValues korrigiert die Imaginaerwerte fuer stetigen Verlauf
% Inputs
% - Eig: Struct mit Vektor Eig.Imag mit zwei Werten der Imaginaerteile
% - buffer struct: Puffer mit letzten Werten fuer Maximum und Minimum
% Outputs
% - Eig: Struct mit angehaengtem, korrigierten Imaginaerteilen
% - buffer struct: Puffer ueberschrieben mit neuen Werten fuer Maximum und Minimum

% Sortiere Imaginaerteil
Eig.ImagSort = sort(Eig.Imag);
% Imaginaeranteil kontinuierlich steigend oder fallend
tmp = Eig.ImagSort(2);
tmpNeg = Eig.ImagSort(1);
% Initialize buffer.Pos/buffer.Neg if they don't exist (safety)
if ~isfield(buffer, 'Pos'), buffer.Pos = 0; end
if ~isfield(buffer, 'Neg'), buffer.Neg = 0; end

if buffer.Pos <= tmp || (abs(tmp) < 10^-5) % Wert uebernehmen
    Eig.ImagCorrected = tmp; % steigender pos. Wert
    Eig.ImagCorrectedNeg = tmpNeg;  % fallender neg. Wert
    buffer.Pos = 0; % Removed based on original logic, these were only placeholders
    buffer.Neg = 0;
else  % 'korrigierten' Wert nehmen fuer kontinuierlichen Verlauf
    Eig.ImagCorrected = 2*buffer.Pos  + tmpNeg; % korrigierter pos. Wert
    Eig.ImagCorrectedNeg  = 2*buffer.Neg + tmp; % korrigierter neg. Wert
end

% Buffer ueberschreiben mit aktuellem Wert
buffer.Pos = max(tmp,buffer.Pos); % Maximum speichern
buffer.Neg = min(tmpNeg,buffer.Neg); % Minimum speichern
end

% 2. DGL System
function dxdpsi = MathieuDGL(psi, x, D, nu_02, nu_C2)
% DGL: phi'' + 2D*phi' + (nu_0^2 + nu_C^2*cos(psi))*phi = 0
% The DGL is converted into a first-order state-space form.
    
    phi = x(1);
    phi_dot = x(2);
    
    % Coefficient of the periodic  stiffness
    K_psi = nu_02 + nu_C2 * cos(psi);
    
    % Second derivative: phi_ddot = -2D*phi_dot - K_psi * phi
    phi_ddot = -2 * D * phi_dot - K_psi * phi;
    
    dxdpsi = [phi_dot; phi_ddot];
end