%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mathieu_Floquet_Verification_m0_argmax.m (Periodicity, Diagnostics & SVG)
%
% Task: Reconstruct x(t) via Peters V(t) with THREE addition factors:
%   (a) Arnold winding number m           (red markers)
%   (b) m = 0, principal exponents        (green markers)
%   (c) argmax-participation factor m_arg (magenta markers, only where
%       it differs from the Arnold value, i.e. cases 2 and 3)
% and validate all against direct ODE integration.
%
% Demonstrates Peters' invariance: the split of a Floquet solution into
% "periodic eigenvector" and "exponential" is arbitrary up to the
% addition factor m, because only the PRODUCT is unique:
%
%   x(t) = V_m(t) * diag(exp(s_m,k * t)) * V0^-1 * x0  =  Phi(t) * x0
%
% with  s_m,1 = sigma + 1i*(omega - m),  V_m(:,1)(t) = V_0(:,1)(t)*e^(+1i*m*t)
%       s_m,2 = sigma + 1i*(omega + m),  V_m(:,2)(t) = V_0(:,2)(t)*e^(-1i*m*t)
%
% The shift e^(+-1i*m*t) cancels between the two factors for ANY m, so
% all reconstructions are identical to each other and to the direct ODE
% solution -- including the argmax-based one, which differs from Arnold
% by a full integer at nu_c^2 = 5 and 8 and is therefore the "wrong"
% frequency label attached to the SAME solution.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;

%% 1. User Settings & Inputs
D = 0.15;                        % Damping ratio
nuInputList = [0.5, 5, 8];       % Amplification factors nu_c^2 (= nu_0^2)
mFactorList = [0.5, 2, 2.5];     % Arnold addition factors for the shift
mFactorArgmax = [0.5, 1.0, 1.5]; % Addition factors implied by dominant
                                 % modal participation (argmax rule);
                                 % case 1 coincides with Arnold
T = 2*pi;                        % Parametric period (T = 2*pi)
t0 = 0;
tGrid = 0:0.1:2*pi;              % 63 distinct timesteps (1x63 vector)
Nz = 2;                          % Number of state dimensions

% --- Setup for Figure Saving ---
fDir = 'figureFolder';           % Target directory for figures
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
    m_arg    = mFactorArgmax(caseIdx);
    doArgmax = abs(m_arg - m_factor) > eps;   % only where the rules differ

    % Coefficients used DIRECTLY as nu_0^2 = nu_c^2, matching the sweep
    % script (where the ODE coefficient equals the x-axis value) and the
    % operating points at which mFactorList/mFactorArgmax were determined.
    % (Previous version squared nuIn here, which evaluated different points.)
    nu02 = nuIn;
    nuC2 = nuIn;

    fprintf('\n=====================================================\n');
    fprintf('Case %d: nu_c^2 = %.1f | Arnold m = %.1f | argmax m = %.1f\n', ...
        caseIdx, nuIn, m_factor, m_arg);
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

    %% Step B: Eigenanalysis & Characteristic Exponent Calculation
    [V0, Mu] = eig(Monodromy);
    multipliers = diag(Mu);

    % Compute Real and Imaginary parts cleanly to prevent asymmetric splitting
    Eig_Real = (1/T) * log(abs(multipliers));
    Eig_Imag = (1/T) * atan2(imag(multipliers), real(multipliers));

    % (a) Shifted exponents according to Arnold addition factor m
    s_R = zeros(Nz, 1);
    s_R(1) = Eig_Real(1) + 1i * (Eig_Imag(1) - m_factor);
    s_R(2) = Eig_Real(2) + 1i * (Eig_Imag(2) + m_factor);

    % (b) Principal exponents: addition factor m = 0 (no shift at all)
    s_P = Eig_Real + 1i * Eig_Imag;

    % (c) Shifted exponents according to argmax-participation factor
    s_A = zeros(Nz, 1);
    s_A(1) = Eig_Real(1) + 1i * (Eig_Imag(1) - m_arg);
    s_A(2) = Eig_Real(2) + 1i * (Eig_Imag(2) + m_arg);

    fprintf('  Shifted exponents s_R (Arnold m = %.1f):\n', m_factor);
    fprintf('    s_R1 = %.4f %+.4fi\n', real(s_R(1)), imag(s_R(1)));
    fprintf('    s_R2 = %.4f %+.4fi\n', real(s_R(2)), imag(s_R(2)));
    fprintf('  Principal exponents s_P (m = 0):\n');
    fprintf('    s_P1 = %.4f %+.4fi\n', real(s_P(1)), imag(s_P(1)));
    fprintf('    s_P2 = %.4f %+.4fi\n', real(s_P(2)), imag(s_P(2)));
    if doArgmax
        fprintf('  Shifted exponents s_A (argmax m = %.1f):\n', m_arg);
        fprintf('    s_A1 = %.4f %+.4fi\n', real(s_A(1)), imag(s_A(1)));
        fprintf('    s_A2 = %.4f %+.4fi\n', real(s_A(2)), imag(s_A(2)));
    else
        fprintf('  argmax factor equals Arnold factor -> no separate s_A case\n');
    end

    %% Step C: Compute the Time-Varying Eigenvector Matrices V(t)
    V_t  = zeros(Nz, Nz, length(tGrid));   % shifted (Arnold m)
    V0_t = zeros(Nz, Nz, length(tGrid));   % principal (m = 0)
    VA_t = zeros(Nz, Nz, length(tGrid));   % shifted (argmax m)
    for j = 1:length(tGrid)
        tj = tGrid(j);
        V_t(:,:,j)  = Phi_t(:,:,j) * V0 * diag([exp(-s_R(1)*tj), exp(-s_R(2)*tj)]);
        V0_t(:,:,j) = Phi_t(:,:,j) * V0 * diag([exp(-s_P(1)*tj), exp(-s_P(2)*tj)]);
        VA_t(:,:,j) = Phi_t(:,:,j) * V0 * diag([exp(-s_A(1)*tj), exp(-s_A(2)*tj)]);
    end

    % Periodicity of the splits:
    % V0_t is ALWAYS T-periodic since exp(s_P*T) equals the multiplier exactly.
    % Shifted splits are T-periodic only for integer m; for half-integer m
    % they are anti-periodic, V(T) = -V(0) (period-doubled branch).
    periodicityError_m  = norm(V_t(:,:,end)  - V_t(:,:,1),  'fro');
    periodicityError_m0 = norm(V0_t(:,:,end) - V0_t(:,:,1), 'fro');
    fprintf('  ||V(T) - V(0)||_F  (Arnold m = %.1f) = %.6e\n', m_factor, periodicityError_m);
    fprintf('  ||V(T) - V(0)||_F  (m = 0)          = %.6e\n', periodicityError_m0);
    if doArgmax
        periodicityError_mA = norm(VA_t(:,:,end) - VA_t(:,:,1), 'fro');
        fprintf('  ||V(T) - V(0)||_F  (argmax m = %.1f) = %.6e\n', m_arg, periodicityError_mA);
    end

    % Relation between the Arnold and principal splits (transfer factor):
    %   V0_t(:,1,j) = V_t(:,1,j) * exp(-1i*m_factor*tGrid(j))
    %   V0_t(:,2,j) = V_t(:,2,j) * exp(+1i*m_factor*tGrid(j))
    transferErr = 0;
    for j = 1:length(tGrid)
        tj = tGrid(j);
        Vcheck = V_t(:,:,j) * diag([exp(-1i*m_factor*tj), exp(+1i*m_factor*tj)]);
        transferErr = max(transferErr, norm(Vcheck - V0_t(:,:,j), 'fro'));
    end
    fprintf('  Max transfer-factor error ||V_m*diag(e^{-imt},e^{+imt}) - V_0|| = %.6e\n', transferErr);

    %% Write out the factorizations x(t) = V(t)*exp(s t)*V0^-1*x0 at a sample time
    jS = 32;  tS = tGrid(jS);          % sample time near t = pi
    Em  = diag([exp(s_R(1)*tS), exp(s_R(2)*tS)]);
    E0  = diag([exp(s_P(1)*tS), exp(s_P(2)*tS)]);
    Pm  = V_t(:,:,jS)  * Em;           % product, Arnold split
    P0  = V0_t(:,:,jS) * E0;           % product, principal split

    fprintf('\n  --- Factorizations at t = %.4f ---\n', tS);
    fprintf('  Arnold split (m = %.1f):   x(t) = V_m(t) * diag(e^{s_R t}) * V0^-1 * x0\n', m_factor);
    printC('    V_m(t)              ', V_t(:,:,jS));
    printC('    diag(e^{s_R t})     ', Em);
    printC('    V_m(t)*diag(e^{s_R t})', Pm);
    fprintf('  Principal split (m = 0):   x(t) = V_0(t) * diag(e^{s_P t}) * V0^-1 * x0\n');
    printC('    V_0(t)              ', V0_t(:,:,jS));
    printC('    diag(e^{s_P t})     ', E0);
    printC('    V_0(t)*diag(e^{s_P t})', P0);
    fprintf('  Max |product difference| = %.6e   (identical: only the product is unique)\n', ...
            max(abs(Pm(:) - P0(:))));
    if doArgmax
        EA = diag([exp(s_A(1)*tS), exp(s_A(2)*tS)]);
        PA = VA_t(:,:,jS) * EA;        % product, argmax split
        fprintf('  Argmax split (m = %.1f):   x(t) = V_A(t) * diag(e^{s_A t}) * V0^-1 * x0\n', m_arg);
        printC('    V_A(t)              ', VA_t(:,:,jS));
        printC('    diag(e^{s_A t})     ', EA);
        printC('    V_A(t)*diag(e^{s_A t})', PA);
        fprintf('  Max |product difference (argmax vs Arnold)| = %.6e\n', ...
                max(abs(PA(:) - Pm(:))));
    end
    fprintf('\n');

    %% Setup Temporary Cell Arrays for Diagnostic Plotting
    X_floquet_cell  = cell(1, 2);      % Arnold-m reconstruction
    X_floquet0_cell = cell(1, 2);      % m = 0 reconstruction
    X_floquetA_cell = cell(1, 2);      % argmax-m reconstruction

    X_direct_cell   = cell(1, 2);

    %% Step D: Trajectory Verification for Initial Conditions
    for initIdx = 1:length(x0_cases)
        x0 = x0_cases{initIdx};

        X_floquet  = zeros(Nz, length(tGrid));
        X_floquet0 = zeros(Nz, length(tGrid));
        X_floquetA = zeros(Nz, length(tGrid));
        V0inv = inv(V0);

        for j = 1:length(tGrid)
            tj = tGrid(j);
            exp_matrix_m  = diag([exp(s_R(1)*tj), exp(s_R(2)*tj)]);
            exp_matrix_m0 = diag([exp(s_P(1)*tj), exp(s_P(2)*tj)]);
            exp_matrix_mA = diag([exp(s_A(1)*tj), exp(s_A(2)*tj)]);
            X_floquet(:,j)  = V_t(:,:,j)  * exp_matrix_m  * V0inv * x0;
            X_floquet0(:,j) = V0_t(:,:,j) * exp_matrix_m0 * V0inv * x0;
            X_floquetA(:,j) = VA_t(:,:,j) * exp_matrix_mA * V0inv * x0;
        end
        X_floquet  = real(X_floquet);
        X_floquet0 = real(X_floquet0);
        X_floquetA = real(X_floquetA);

        % Direct ground-truth solution via ODE solver
        solDirect = ode45(@(t,x) MathieuDGL_task(t, x, D, nu02, nuC2), [t0, T], x0, opts);
        X_direct = deval(solDirect, tGrid);

        X_floquet_cell{initIdx}  = X_floquet;
        X_floquet0_cell{initIdx} = X_floquet0;
        X_floquetA_cell{initIdx} = X_floquetA;
        X_direct_cell{initIdx}   = X_direct;

        % Calculate global trajectory errors
        maxErr_m  = max(vecnorm(X_floquet  - X_direct,   2, 1));
        maxErr_m0 = max(vecnorm(X_floquet0 - X_direct,   2, 1));
        maxErr_mA = max(vecnorm(X_floquetA - X_direct,   2, 1));
        fprintf('  IC [%d; %d] -> Max Err (Arnold m = %.1f) = %.4e | (m = 0) = %.4e | (argmax m = %.1f) = %.4e\n', ...
                x0(1), x0(2), m_factor, maxErr_m, maxErr_m0, m_arg, maxErr_mA);

        %% Step E: Save to global struct matrix
        AllResults(structIdx).nu = nuIn;
        AllResults(structIdx).m = m_factor;
        AllResults(structIdx).m_arg = m_arg;
        AllResults(structIdx).x0 = x0;
        AllResults(structIdx).s_R = s_R;
        AllResults(structIdx).s_P = s_P;
        AllResults(structIdx).s_A = s_A;
        AllResults(structIdx).V_t = V_t;
        AllResults(structIdx).V0_t = V0_t;
        AllResults(structIdx).VA_t = VA_t;
        AllResults(structIdx).X_floquet = X_floquet;
        AllResults(structIdx).X_floquet0 = X_floquet0;
        AllResults(structIdx).X_floquetA = X_floquetA;
        AllResults(structIdx).X_direct = X_direct;
        AllResults(structIdx).maxErr = maxErr_m;
        AllResults(structIdx).maxErr_m0 = maxErr_m0;
        AllResults(structIdx).maxErr_mA = maxErr_mA;

        structIdx = structIdx + 1;
    end

    %% Step F: Diagnostic Visualization Window
    hFig = figure('Name', sprintf('Case %d Validation Plot', caseIdx), 'Color', [1 1 1]);



tlo = tiledlayout(2,1,'TileSpacing','tight','Padding','compact');

    % --- Subplot 1: Displacement ---
    nexttile;
    hold on; grid on;
    h1 = plot(tGrid, X_direct_cell{1}(1,:), 'b-', 'LineWidth', 2.0);
    % m = 0 reconstruction: green markers, drawn first (largest, underneath)
    h5 = plot(tGrid, X_floquet0_cell{1}(1,:), 'gs', 'MarkerSize', 10, 'LineWidth', 1.0);
    if doArgmax
        % argmax reconstruction: magenta markers, in between
        h7 = plot(tGrid, X_floquetA_cell{1}(1,:), 'm^', 'MarkerSize', 7.5, 'LineWidth', 1.0);
    end
    h2 = plot(tGrid, X_floquet_cell{1}(1,:),  'ro', 'MarkerSize', 5,  'LineWidth', 1.2);
    h3 = plot(tGrid, X_direct_cell{2}(1,:), 'b--', 'LineWidth', 1.5);
    h6 = plot(tGrid, X_floquet0_cell{2}(1,:), 'gd', 'MarkerSize', 10, 'LineWidth', 1.0);
    if doArgmax
        h8 = plot(tGrid, X_floquetA_cell{2}(1,:), 'mv', 'MarkerSize', 7.5, 'LineWidth', 1.0);
    end
    h4 = plot(tGrid, X_floquet_cell{2}(1,:),  'rx', 'MarkerSize', 6,  'LineWidth', 1.2);

    ylabel('Displacement \phi(t)', 'FontSize', 11, 'FontWeight', 'bold');
    if doArgmax
        title(sprintf('Mathieu \\nu_c^2 = %.1f: red: m = %.1f, green: m = 0, magenta: m = %.1f', ...
              nuIn, m_factor, m_arg), 'FontSize', 12, 'FontWeight', 'bold');
        legend([h1, h2, h5, h7, h3, h4, h6, h8], ...
               'Direct ODE [1;0]', sprintf('Floquet m=%.1f [1;0]', m_factor), ...
               'Floquet m=0 [1;0]', sprintf('Floquet m=%.1f [1;0]', m_arg), ...
               'Direct ODE [0;1]', sprintf('Floquet m=%.1f [0;1]', m_factor), ...
               'Floquet m=0 [0;1]', sprintf('Floquet m=%.1f [0;1]', m_arg), ...
               'Location', 'best', 'NumColumns', 2);
    else
        title(sprintf('Mathieu Diagnostics: \\nu_c^2 = %.1f | red: m = %.1f, green: m = 0 (argmax = Arnold)', ...
              nuIn, m_factor), 'FontSize', 12, 'FontWeight', 'bold');
        legend([h1, h2, h5, h3, h4, h6], ...
               'Direct ODE [1;0]', sprintf('Floquet m=%.1f [1;0]', m_factor), 'Floquet m=0 [1;0]', ...
               'Direct ODE [0;1]', sprintf('Floquet m=%.1f [0;1]', m_factor), 'Floquet m=0 [0;1]', ...
               'Location', 'best');
    end
    set(gca, 'XTick', 0:pi/2:2*pi, 'XTickLabel', {});
    xlim([0 2*pi]);

    % --- Subplot 2: Velocity ---
    nexttile;
    hold on; grid on;
    plot(tGrid, X_direct_cell{1}(2,:), 'b-', 'LineWidth', 2.0);
    plot(tGrid, X_floquet0_cell{1}(2,:), 'gs', 'MarkerSize', 10, 'LineWidth', 1.0);
    if doArgmax
        plot(tGrid, X_floquetA_cell{1}(2,:), 'm^', 'MarkerSize', 7.5, 'LineWidth', 1.0);
    end
    plot(tGrid, X_floquet_cell{1}(2,:),  'ro', 'MarkerSize', 5,  'LineWidth', 1.2);
    plot(tGrid, X_direct_cell{2}(2,:), 'b--', 'LineWidth', 1.5);
    plot(tGrid, X_floquet0_cell{2}(2,:), 'gd', 'MarkerSize', 10, 'LineWidth', 1.0);
    if doArgmax
        plot(tGrid, X_floquetA_cell{2}(2,:), 'mv', 'MarkerSize', 7.5, 'LineWidth', 1.0);
    end
    plot(tGrid, X_floquet_cell{2}(2,:),  'rx', 'MarkerSize', 6,  'LineWidth', 1.2);

    xlabel('Time t [rad]', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('Velocity d\phi/dt', 'FontSize', 11, 'FontWeight', 'bold');
    set(gca, 'XTick', 0:pi/2:2*pi, ...
             'XTickLabel', {'0', '\pi/2', '\pi', '3\pi/2', '2\pi'});
    xlim([0 2*pi]);

    %% Save Figure directly to figureFolder as an SVG
    % (dots inside the name replaced BEFORE appending the extension)
    baseName = sprintf('Mathieu_Verification_m0_argmax_Case_%d_nu_%.1f_m_%.1f_mA_%.1f', ...
        caseIdx, nuIn, m_factor, m_arg);
    baseName = strrep(baseName, '.', 'dot');
    svgFileFull = fullfile(fDir, [baseName, '.svg']);

    print(hFig, svgFileFull, '-dsvg'); % Export as vector graphic file
end

% Save verified results matrix
save('Mathieu_Floquet_Verification_m0_argmax_Results.mat', 'AllResults');
fprintf('\n=====================================================\n');
disp('Verification completed: m = 0 (green), Arnold m (red) and argmax m');
disp('(magenta, cases 2 & 3) yield identical trajectories, since only');
disp('the product V(t)*e^{st} is unique.');

%% Local Function: First-Order Damped Mathieu State Space System
function dx = MathieuDGL_task(t, x, D, nu02, nuC2)
    dx = zeros(2,1);
    dx(1) = x(2);
    dx(2) = -(nu02 + nuC2*cos(t))*x(1) - 2*D*x(2);
end

%% Local Function: Pretty-print a 2x2 complex matrix on two lines
function printC(label, M)
    fprintf('%s = [% .4f%+.4fi   % .4f%+.4fi ;\n', label, ...
            real(M(1,1)), imag(M(1,1)), real(M(1,2)), imag(M(1,2)));
    fprintf('%s    % .4f%+.4fi   % .4f%+.4fi ]\n', blanks(length(label)), ...
            real(M(2,1)), imag(M(2,1)), real(M(2,2)), imag(M(2,2)));
end