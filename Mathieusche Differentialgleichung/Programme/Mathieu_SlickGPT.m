%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation of the Strutt Diagrams within the bounds of nu_02 and nu_C2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;

% Configuration parameters
loadMat = 0;        % Load mat-file if results with same D already exist
SW = 0.5;           % Step width
unt0 = 0;          % Lower bound for nu_02
untC = 0;          % Lower bound for nu_C2
ob0 = 9;           % Upper bound for nu_02
obC = 9;           % Upper bound for nu_C2
fDir = 'figureFolder1';   % Folder for figures
dDir = 'dataFolder';       % Folder for mat-files
excelDir = 'dataFolder1';  % Folder for exporting Excel files
Nz = 2;                     % Number of equations in the ODE system
DVec = 0.15;                % D parameter

% Create necessary directories
if ~isfolder(fDir), mkdir(fDir); end
if ~isfolder(dDir), mkdir(dDir); end
if ~isfolder(excelDir), mkdir(excelDir); end

t0 = 0.0;                   % Initial time
T = 2 * pi;                 % Final time

for dIdx = 1:length(DVec)
    D = DVec(dIdx);
    
    % Initial conditions
    Diagonal = eye(Nz);      % Identity matrix for initial conditions
    matName = sprintf('STRUTTscheKarte_D%2.1e_SW%2.1e.mat', D, SW);
    
    % Generate mat-file name based on the conditions
    if unt0 == 0
        matName = strrep(matName, '.mat', '_unt0.mat');
    end
    fileName = fullfile(dDir, matName);
    
    % Preparation of arrays
    lenNu02 = length(unt0:SW:ob0);
    lenNuC2 = length(untC:SW:obC);
    lenNu = lenNu02 * lenNuC2;
    lenNuDiag = min(lenNu02, lenNuC2);
    
    % Initialize storage matrices
    Monodromie = zeros(Nz);
    CharEx = zeros(lenNuDiag, Nz * 4 + 2);
    nAddVector = nan(lenNuDiag, 1);
    plotwertstabil = zeros(lenNu, 3);
    
    nuCSwitchVec = [0.2, 1, 1.5^2, 2^2, 2.5^2] - 0.15;
    n = 0;
    cntN = 1;
    buffer = struct('Pos', 0, 'Neg', 0);
    
    % Load or compute data
    if exist(fileName, 'file') == 2 && loadMat == 1
        load(fileName, 'Char*', 'plotwert*', 'nu*', 'nAddVector');
    else
        % Variable indices initialization for characteristic multipliers
        lidx = 1; % Index for stable combinations
        oidx = 1; % Index for characteristic exponents
        
        % Loop through nu_02 and nu_C2 combinations
        for nu_02 = unt0:SW:ob0
            for nu_C2 = untC:SW:obC
                options = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
                for k = 1:Nz
                    % Solve ODEs and store the solution in Monodromy matrix
                    sol = ode45(@(psi, x) MathieuDGL(psi, x, D, nu_02, nu_C2), [t0, T], Diagonal(:, k), options);
                    MonoVek = deval(sol, T);
                    Monodromie(:, k) = MonoVek;
                end
                
                % Calculate characteristic multipliers (eigenvalues of Monodromy matrix)
                eP = eig(Monodromie);
                
                % Process characteristic exponents
                if nu_C2 == nu_02
                    Eig = processEigenvalues(eP);
                    [Eig, buffer] = correctImagValues(Eig, buffer);
                    
                    nAdd = updateNValue(Eig, nu_C2, nuCSwitchVec, cntN, T);

                    ImagEigSortN = [Eig.ImagCorrectedNeg, Eig.ImagCorrected] + [-nAdd,nAdd];
                    CharEx(oidx,:) = [nu_02, nu_C2, Eig.Real', min(Eig.ImagSort), max(Eig.ImagSort), ImagEigSortN, eP'];
                 
                    nAddVector(oidx) = nAdd;
                    oidx = oidx + 1;
                end
                
                % Record stable combinations
                if max(abs(eP)) < 1
                    plotwertstabil(lidx, :) = [nu_02, nu_C2, 1]; % 1 indicates stability
                end
                lidx = lidx + 1;
            end
        end
        
        % Remove empty rows and save results
        plotwertstabil = plotwertstabil(any(plotwertstabil, 2), :);
        save(fileName, 'Char*', 'plotwert*', 'nu*', 'nAddVector');
    end

    %% Export results to Excel
    try
        exportToExcel(fileName, CharEx);
    catch
        warning('Failed to write Excel file: %s', fileName);
    end

    %% Graphical representation
    plotResults(plotwertstabil, CharEx, dIdx, matName, fDir, D);
end

function Eig = processEigenvalues(eP)
    % Calculate the characteristic eigenvalues
    Eig.Real = log(abs(eP)) / (2 * pi);
    Eig.Imag = angle(eP) / (2 * pi);
    Eig.ImagSort = sort(Eig.Imag);
end

function nAdd = updateNValue(Eig, nu_C2, nuCSwitchVec, cntN, T)
    % Update the n value based on switches in the characteristic exponents
    nAdd = 0;
    if cntN <= length(nuCSwitchVec) && nu_C2 > nuCSwitchVec(cntN)
        if abs(Eig.Real(1) - Eig.Real(2)) > eps
            nAdd = nAdd + 0.5; % Increment only on switch
            cntN = cntN + 1;
        end
    end
    nAdd = nAdd * 2 * pi / T;
end

function exportToExcel(fileName, CharEx)
    % Function to export computed characteristics to an Excel file
    CharExTable = array2table(CharEx(:, [1:4, 7, 8]), 'VariableNames', ...
                               {'nu02', 'nu_C2', 'Eig.Real1', 'Eig.Real2', ...
                               'ImagEig1', 'ImagEig2'});
    charExFileName = strrep(fileName, '.mat', 'CharExAll.xlsx');
    writetable(CharExTable, charExFileName);
end

function plotResults(plotwertstabil, CharEx, dIdx, matName, fDir, D)
    % Function to create plots of the results
    fs = 10.5; % Fontsize
    hf = figure(dIdx);
    hf.Position = [10, 10, 600, 600];
    
    try tiledlayout(11, 1); catch, end % tiledlayout not available before R2007
    h(1) = nexttile([7 1]); % Strutt diagram
    scatter(plotwertstabil(:, 1), plotwertstabil(:, 2), 5, 'filled');
    ylabel('$\rm{Parameter} \; \nu_C^2 \;\;\rm{[-]}$', 'interpreter', 'latex', 'Position', [-0.5 4.5], 'FontSize', fs + 2);
    grid on;

    % Real parts of characteristic exponents
    h(2) = nexttile([2 1]);
    plot(CharEx(:, 1), CharEx(:, 3:4));
    ylabel('$Re(s_R) \;\; \rm{[-]}$', 'interpreter', 'latex', 'FontSize', fs);
    grid on;

    % Imaginary parts of characteristic exponents
    h(3) = nexttile([2 1]);
    plot(CharEx(:, 1), CharEx(:, 7:8)); hold on;
    ylim(gca, [-0.25, 3]);
    xlabel('$\rm{Parameter} \; \nu_0^2 \;\;\rm{[-]}$', 'interpreter', 'latex', 'FontSize', fs + 2);
    ylabel('$Im(s_R) \;\; \rm{[-]}$', 'interpreter', 'latex', 'FontSize', fs);
    linkaxes(h, 'x');

    for idxH = 1:length(h)
        set(h(idxH), 'TickLabelInterpreter', 'Latex', 'FontSize', fs);
    end

    % Save figure as PNG
    pngname = fullfile(fDir, strrep(matName, '.mat', ''));
    print(pngname, '-dpng');

    % Second figure for detailed view
    hf1 = figure(dIdx + 100);
    h(1) = tiledlayout(2, 1);
    nexttile;
    plot(CharEx(:, 1), CharEx(:, 3:4)); grid on;
    title({['$\ddot{\phi} + 2D \dot{\phi} + (\nu^2_0 + \nu^2_C \cos(\psi))\phi = 0$, D = ', num2str(D, 2), ...
            ', $\nu_C = \nu_0 $'] ...
           '$\Re(s_{Ri}) = 1/2\pi\ln(|\mu_{Ri}|)$, $\mu_{Ri}$: Eigenvalue monodromy matrix'}, ...
           'interpreter', 'latex', 'FontSize', fs + 2);
    ylabel('$\Re(s_R) \;\; \rm{[-]}$', 'interpreter', 'latex', 'FontSize', fs + 2);

    nexttile;
    plot(CharEx(:, 1), CharEx(:, 7:8)); hold on; grid on;
    xlabel('$\rm{Parameter} \; \nu_C^2 \;\;\rm{[-]}$', 'interpreter', 'latex', 'FontSize', fs + 2);
    ylabel('$\Im(s_R) \;\; \rm{[-]}$', 'interpreter', 'latex', 'FontSize', fs);
end

function [Eig, buffer] = correctImagValues(Eig, buffer)
    % correctImagValues corrects the imaginary values for continuous behavior
    % Inputs
    % - Eig: Struct with a vector Eig.Imag containing two imaginary parts
    % - buffer: Struct with last values for max and min
    % Outputs
    % - Eig: Struct with connected corrected imaginary parts
    % - buffer: Updated buffer with new max and min values

    % Check input format
    if nargin~= 2
        error('Two inputs are expected.');
    end
    if max(size(Eig.Imag)) ~= 2 || min(size(Eig.Imag)) ~= 1
        error('The current imaginary parts of an eigenvalue pair are expected.');
    end

    % Sort imaginary part
    Eig.ImagSort = sort(Eig.Imag);
    tmp = Eig.ImagSort(2);
    tmpNeg = Eig.ImagSort(1);

    % Correct imaginary values based on the buffer
    if buffer.Pos <= tmp || (abs(tmp) < 1e-5) % Adopt the value
        Eig.ImagCorrected = tmp; % increasing positive value
        Eig.ImagCorrectedNeg = tmpNeg; % decreasing negative value
        buffer.Pos = 0;
        buffer.Neg = 0;
    else % Take the corrected value for continuous behavior
        Eig.ImagCorrected = 2 * buffer.Pos + tmpNeg; % corrected pos. value
        Eig.ImagCorrectedNeg = 2 * buffer.Neg + tmp; % corrected neg. value
    end

    % Update the buffer with the latest values
    buffer.Pos = max(tmp, buffer.Pos); % Save maximum
    buffer.Neg = min(tmpNeg, buffer.Neg); % Save minimum
end
