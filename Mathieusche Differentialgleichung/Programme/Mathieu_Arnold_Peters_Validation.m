%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mathieu_Floquet_Verification.m (Fixed Periodicity, Diagnostics & SVG Save)
%
% Task: Reconstruct x(t) (2x63) via Peters V(t) and Arnold s_R, 
% then validate against direct ODE integration. Saves plots as SVG.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;

%% 1. User Settings & Inputs
D = 0.15;                       % Damping ratio
nuInputList = [0.5, 5, 8];      % Frequency parameters (nu_0 = nu_C)
mFactorList = [0.5, 2, 2.5];    % Arnold addition factors for the shift
T = 2*pi;                       % Parametric period (T = 2*pi)
t0 = 0;
tGrid = 0:0.1:2*pi;             % 63 distinct timesteps (1x63 vector)
Nz = 2;                         % Number of state dimensions

% --- Setup for Figure Saving ---
fDir = 'figureFolder';          % Target directory for figures
if ~exist(fDir, 'dir')
    mkdir(fDir);
end

% Complete list of required initial conditions
x0_cases = {[1; 0], [0; 1]};
opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);

AllResults = struct();
structIdx = 1;

%% 2. Execution Loop across Parameter Cases
for caseIdx = 1:length(nuInputList)
    nuIn = nuInputList(caseIdx);
    m_factor = mFactorList(caseIdx);
    
    nu02 = nuIn^2;
    nuC2 = nuIn^2;
    
    fprintf('\n=====================================================\n');
    fprintf('Case %d: nu_0 = nu_C = %.1f | Arnold m = %.1f\n', caseIdx, nuIn, m_factor);
    fprintf('=====================================================\n');
    
    %% Step A: Compute the fundamental transition matrix Phi(t)
    Phi_t = zeros(Nz, Nz, length(tGrid));
    Monodromy = zeros(Nz, Nz);
    I2 = eye(Nz);
    
    for k = 1:Nz
        % Integrate standard basis unit vectors over the full period
        solBasis = ode45(@(t,x) MathieuDGL_task(t, x, D, nu02, nuC2), [t0, T], I2(:,k), opts);
        for j = 1:length(tGrid)
            Phi_t(:,k,j) = deval(solBasis, tGrid(j));
        end
        Monodromy(:,k) = deval(solBasis, T); 
    end
    
    %% Step B: Eigenanalysis & Arnold Exponent Matrix Calculation
    [V0, Mu] = eig(Monodromy);
    multipliers = diag(Mu);
    
    % Compute Real and Imaginary parts cleanly to prevent asymmetric splitting
    Eig_Real = (1/T) * log(abs(multipliers));
    Eig_Imag = (1/T) * atan2(imag(multipliers), real(multipliers));
    
    % Formulate the characteristic exponents according to Arnold shifts
    s_R = zeros(Nz, 1);
    s_R(1) = Eig_Real(1) + 1i * (Eig_Imag(1) - m_factor);
    s_R(2) = Eig_Real(2) + 1i * (Eig_Imag(2) + m_factor);
    
    fprintf('  Arnold Exponents s_R:\n');
    fprintf('    s_R1 = %.4f + %.4fi\n', real(s_R(1)), imag(s_R(1)));
    fprintf('    s_R2 = %.4f + %.4fi\n', real(s_R(2)), imag(s_R(2)));
    
    %% Step C: Compute Time-Varying Eigenvector Matrix V(t)
    V_t = zeros(Nz, Nz, length(tGrid));
    for j = 1:length(tGrid)
        tj = tGrid(j);
        V_t(:,:,j) = Phi_t(:,:,j) * V0 * diag([exp(-s_R(1)*tj), exp(-s_R(2)*tj)]);
    end
    
    % Verify structural periodicity of Peters' matrix
    periodicityError = norm(V_t(:,:,end) - V_t(:,:,1), 'fro');
    fprintf('  ||V(T) - V(0)||_F Periodicity Error = %.6e\n', periodicityError);
    
    %% Setup Temporary Cell Arrays for Diagnostic Plotting
    X_floquet_cell = cell(1, 2);
    X_direct_cell  = cell(1, 2);
    
    %% Step D: Trajectory Verification for Initial Conditions
    for initIdx = 1:length(x0_cases)
        x0 = x0_cases{initIdx};
        
        X_floquet = zeros(Nz, length(tGrid));
        V0inv = inv(V0);
        
        for j = 1:length(tGrid)
            tj = tGrid(j);
            exp_matrix = diag([exp(s_R(1)*tj), exp(s_R(2)*tj)]);
            X_floquet(:,j) = V_t(:,:,j) * exp_matrix * V0inv * x0;
        end
        X_floquet = real(X_floquet); 
        
        % Direct ground-truth solution via ODE solver
        solDirect = ode45(@(t,x) MathieuDGL_task(t, x, D, nu02, nuC2), [t0, T], x0, opts);
        X_direct = deval(solDirect, tGrid);
        
        X_floquet_cell{initIdx} = X_floquet;
        X_direct_cell{initIdx}  = X_direct;
        
        % Calculate global trajectory error
        X_error = X_floquet - X_direct;
        maxErr = max(vecnorm(X_error, 2, 1));
        fprintf('  Initial Condition [%d; %d] -> Max Abs Error = %.6e\n', x0(1), x0(2), maxErr);
        
        %% Step E: Save to global struct matrix
        AllResults(structIdx).nu = nuIn;
        AllResults(structIdx).m = m_factor;
        AllResults(structIdx).x0 = x0;
        AllResults(structIdx).s_R = s_R;
        AllResults(structIdx).V_t = V_t;
        AllResults(structIdx).X_floquet = X_floquet;
        AllResults(structIdx).X_direct = X_direct;
        AllResults(structIdx).maxErr = maxErr;
        
        structIdx = structIdx + 1;
    end
    
    %% Step F: Diagnostic Visualization Window
    hFig = figure('Name', sprintf('Case %d Validation Plot', caseIdx), 'Color', [1 1 1]);
    
    % --- Subplot 1: Displacement ---
    subplot(2,1,1);
    hold on; grid on;
    h1 = plot(tGrid, X_direct_cell{1}(1,:), 'b-', 'LineWidth', 2.0);
    h2 = plot(tGrid, X_floquet_cell{1}(1,:), 'ro', 'MarkerSize', 5, 'LineWidth', 1.2);
    h3 = plot(tGrid, X_direct_cell{2}(1,:), 'b--', 'LineWidth', 1.5);
    h4 = plot(tGrid, X_floquet_cell{2}(1,:), 'rx', 'MarkerSize', 6, 'LineWidth', 1.2);
    
    ylabel('Displacement \phi(t)', 'FontSize', 11, 'FontWeight', 'bold');
    title(sprintf('Mathieu Diagnostics: \\nu_C = %.1f, Addition Factor m = %.1f', nuIn, m_factor), ...
          'FontSize', 12, 'FontWeight', 'bold');
    legend([h1, h2, h3, h4], ...
           'Direct ODE [1;0]', 'Floquet Recon [1;0]', ...
           'Direct ODE [0;1]', 'Floquet Recon [0;1]', ...
           'Location', 'best');
    set(gca, 'XTick', 0:pi/2:2*pi, 'XTickLabel', {}); 
    xlim([0 2*pi]);
    
    % --- Subplot 2: Velocity ---
    subplot(2,1,2);
    hold on; grid on;
    plot(tGrid, X_direct_cell{1}(2,:), 'b-', 'LineWidth', 2.0);
    plot(tGrid, X_floquet_cell{1}(2,:), 'ro', 'MarkerSize', 5, 'LineWidth', 1.2);
    plot(tGrid, X_direct_cell{2}(2,:), 'b--', 'LineWidth', 1.5);
    plot(tGrid, X_floquet_cell{2}(2,:), 'rx', 'MarkerSize', 6, 'LineWidth', 1.2);
    
    xlabel('Time t [rad]', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('Velocity d\phi/dt', 'FontSize', 11, 'FontWeight', 'bold');
    set(gca, 'XTick', 0:pi/2:2*pi, ...
             'XTickLabel', {'0', '\pi/2', '\pi', '3\pi/2', '2\pi'});
    xlim([0 2*pi]);
    
    %% Save Figure directly to figureFolder as an SVG
    svgName = sprintf('Mathieu_Verification_Case_%d_nu_%.1f_m_%.1f.svg', caseIdx, nuIn, m_factor);
    svgName = strrep(svgName, '.', 'dot'); % Prevent extension confusion from dots in floats
    svgName = strrep(svgName, 'dotyvgz', '.svg'); % Re-establish the pure file extension block
    svgFileFull = fullfile(fDir, svgName);
    
    print(hFig, svgFileFull, '-dsvg'); % Export as vector graphic file
end

% Save verified results matrix
save('Mathieu_Floquet_Verification_Results.mat', 'AllResults');
fprintf('\n=====================================================\n');
disp('Verification completed successfully! Plots saved as SVG in figureFolder.');

%% Local Function: First-Order Damped Mathieu State Space System
function dx = MathieuDGL_task(t, x, D, nu02, nuC2)
    dx = zeros(2,1);
    dx(1) = x(2); 
    dx(2) = -(nu02 + nuC2*cos(t))*x(1) - 2*D*x(2); 
end