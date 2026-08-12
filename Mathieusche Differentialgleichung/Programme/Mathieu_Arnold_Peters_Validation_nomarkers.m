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
Tvec = [2*pi, 4*pi];                        % Parametric period (T = T)
T = 2*pi;
t0 = 0;

Nz = 2;                          % Number of state dimensions

% --- Setup for Figure Saving ---
fDir = 'figureFolder1';           % Target directory for figures
if ~exist(fDir, 'dir')
    mkdir(fDir);
end

% Complete list of required initial conditions
x0_cases = {[1; 0], [0; 1]};
opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);

AllResults = struct();
structIdx = 1;

%% 2. Execution Loop across Parameter Cases
for TvecIdx = 2 %1: length(Tvec)
    Tlen = Tvec(TvecIdx);
    tGrid = 0:0.1:Tlen;              % 63 distinct timesteps (1x63 vector)
    
    % Generate correct x-axis labels based on period Tlen
    ticks = 0:pi/2:Tlen;
    strX = cell(1, length(ticks));
    for i = 1:length(ticks)
        k = round(ticks(i) / (pi/2));
        if k == 0
            strX{i} = '0';
        elseif mod(k, 2) == 0
            strX{i} = sprintf('%d\\pi', k/2);
        else
            strX{i} = sprintf('%d\\pi/2', k);
        end
    end
    
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
            solBasis = ode45(@(t,x) MathieuDGL_task(t, x, D, nu02, nuC2), [t0, Tlen], I2(:,k), opts);
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
            solDirect = ode45(@(t,x) MathieuDGL_task(t, x, D, nu02, nuC2), [t0, Tlen], x0, opts);
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

        % Grayscale-safe styling: layered widths, no markers.
        % Direct ODE: widest, lightest, underneath. Arnold m: thinnest, black, on top.
        colDirect = [0.72 0.72 0.72];   % light gray
        colM0     = [0.45 0.45 0.45];   % medium gray
        colArg    = [0.20 0.20 0.20];   % dark gray
        colArnold = [0    0    0   ];   % black

        % Grayscale-safe colors with distinct hues (viridis, yellow end excluded)
        colDirect = [0.478 0.821 0.318];   % light green  (viridis ~0.85, brightest)
        colM0     = [0.128 0.647 0.523];   % teal         (viridis ~0.60)
        colArg    = [0.268 0.409 0.559];   % blue         (viridis ~0.35)
        colArnold = [0.267 0.049 0.335];   % dark purple  (viridis ~0.02)


        cl = parula(256);
        colDirect = cl(200,:);   % green-cyan, brightest used
        colM0     = cl(140,:);   % cyan-blue
        colArg    = cl( 70,:);   % blue
        colArnold = cl(  1,:);   % dark blue
        lwDirect = 4;  lwM0 = 2.2;  lwArg = 1.4;  lwArnold = 2;

        tlo = tiledlayout(2,1,'TileSpacing','tight','Padding','compact');

        % --- Subplot 1: Displacement ---
        nexttile;
        hold on; grid on;
        h1 = plot(tGrid, X_direct_cell{1}(1,:),   '-',  'LineWidth', lwDirect, 'Color', colDirect);
        h5 = plot(tGrid, X_floquet0_cell{1}(1,:), '--', 'LineWidth', lwM0,     'Color', colM0);
        if doArgmax
            h7 = plot(tGrid, X_floquetA_cell{1}(1,:), '-.', 'LineWidth', lwArg, 'Color', colArg);
        end
        h2 = plot(tGrid, X_floquet_cell{1}(1,:),  ':',  'LineWidth', lwArnold, 'Color', colArnold);

        h3 = plot(tGrid, X_direct_cell{2}(1,:),   '-',  'LineWidth', lwDirect, 'Color', colDirect);
        h6 = plot(tGrid, X_floquet0_cell{2}(1,:), '--', 'LineWidth', lwM0,     'Color', colM0);
        if doArgmax
            h8 = plot(tGrid, X_floquetA_cell{2}(1,:), '-.', 'LineWidth', lwArg, 'Color', colArg);
        end
        h4 = plot(tGrid, X_floquet_cell{2}(1,:),  ':',  'LineWidth', lwArnold, 'Color', colArnold);

        ylabel('Displacement \phi(t)', 'FontSize', 11, 'FontWeight', 'bold');
        if doArgmax
            title(sprintf('Mathieu \\nu_c^2 = %.1f: Arnold m = %.1f, argmax m = %.1f, m = 0', ...
                nuIn, m_factor, m_arg), 'FontSize', 12, 'FontWeight', 'bold');
            legend([h1, h5, h7, h2], ...
                'Direct ODE', 'Floquet m=0', ...
                sprintf('Floquet m=%.1f (argmax)', m_arg), ...
                sprintf('Floquet m=%.1f (Arnold)', m_factor), ...
                'Location', 'best');
        else
            title(sprintf('Mathieu Diagnostics: \\nu_c^2 = %.1f | m = %.1f (argmax = Arnold)', ...
                nuIn, m_factor), 'FontSize', 12, 'FontWeight', 'bold');
            legend([h1, h5, h2], ...
                'Direct ODE', 'Floquet m=0', sprintf('Floquet m=%.1f', m_factor), ...
                'Location', 'best');
        end
        set(gca, 'XTick', 0:pi/2:Tlen, 'XTickLabel', {});
        xlim([0 Tlen]);

        % --- Subplot 2: Velocity ---
        nexttile;
        hold on; grid on;
        plot(tGrid, X_direct_cell{1}(2,:),   '-',  'LineWidth', lwDirect, 'Color', colDirect);
        plot(tGrid, X_floquet0_cell{1}(2,:), '--', 'LineWidth', lwM0,     'Color', colM0);
        if doArgmax
            plot(tGrid, X_floquetA_cell{1}(2,:), '-.', 'LineWidth', lwArg, 'Color', colArg);
        end
        plot(tGrid, X_floquet_cell{1}(2,:),  ':',  'LineWidth', lwArnold, 'Color', colArnold);

        plot(tGrid, X_direct_cell{2}(2,:),   '-',  'LineWidth', lwDirect, 'Color', colDirect);
        plot(tGrid, X_floquet0_cell{2}(2,:), '--', 'LineWidth', lwM0,     'Color', colM0);
        if doArgmax
            plot(tGrid, X_floquetA_cell{2}(2,:), '-.', 'LineWidth', lwArg, 'Color', colArg);
        end
        plot(tGrid, X_floquet_cell{2}(2,:),  ':',  'LineWidth', lwArnold, 'Color', colArnold);

        xlabel('Time t [rad]', 'FontSize', 11, 'FontWeight', 'bold');
        ylabel('Velocity d\phi/dt', 'FontSize', 11, 'FontWeight', 'bold');
        set(gca, 'XTick', 0:pi/2:Tlen, ...
            'XTickLabel', strX);
        xlim([0 Tlen]);


        %% Save Figure directly to figureFolder as an SVG
        % (dots inside the name replaced BEFORE appending the extension)
        baseName = sprintf('Mathieu_Verif_NoMarkers_Case_%d_nu_%.1f_m_%.1f_mA_%.1f', ...
            caseIdx, nuIn, m_factor, m_arg);
        baseName = strrep(baseName, '.', 'dot');
        svgFileFull = fullfile(fDir, [baseName, '.svg']);

        print(hFig, svgFileFull, '-dsvg'); % Export as vector graphic file



        %% Step G: Time evolution of V(t) and w(t) = diag(e^{s t}) * V0^-1 * x0
        % x(t) = V_m(t) * w_m(t) with w_m(t) = diag(e^{s_m,k t}) * V0^-1 * x0.
        % The e^{+-i m t} shift migrates between V and w; each factor changes
        % with m, only the product does not.
        icIdx = 1;                          % initial condition to visualize
        x0G = x0_cases{icIdx};
        c0 = V0 \ x0G;                      % constant modal coefficients V0^-1*x0

        nT = length(tGrid);
        wR = zeros(2, nT); w0 = zeros(2, nT); wA = zeros(2, nT);
        for j = 1:nT
            tj = tGrid(j);
            wR(:,j) = [exp(s_R(1)*tj); exp(s_R(2)*tj)] .* c0;
            w0(:,j) = [exp(s_P(1)*tj); exp(s_P(2)*tj)] .* c0;
            wA(:,j) = [exp(s_A(1)*tj); exp(s_A(2)*tj)] .* c0;
        end

        % Colors per m-variant (same scheme as Step F); column 2 lighter
        colVar = {colM0, colArnold, colArg};              % m=0, Arnold, argmax
        varName = {sprintf('m=0'), sprintf('m=%.1f (Arnold)', m_factor), ...
            sprintf('m=%.1f (argmax)', m_arg)};
        Vset = {V0_t, V_t, VA_t};
        wset = {w0, wR, wA};
        nVar = 2 + doArgmax;                              % skip argmax if identical
        lighten = @(c) 1 - 0.55*(1 - c);                  % col-2 shade

        hFigG = figure('Name', sprintf('Case %d Factor Evolution', caseIdx), ...
            'Color', [1 1 1], 'Position', [100 100 900 700]);
        tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

        % --- Tiles 1 & 2: rows of V(t) (position row, velocity row) ---
        rowLabel = {'V(1,:) position row', 'V(2,:) velocity row'};
        for row = 1:2
            nexttile; hold on; grid on;
            hLeg = gobjects(1, nVar);
            for v = 1:nVar
                Vv = Vset{v};
                % column 1: full color, column 2: lighter/thinner
                hLeg(v) = plot(tGrid, real(squeeze(Vv(row,1,:))), '-',  ...
                    'LineWidth', 1.6, 'Color', colVar{v});
                plot(tGrid, imag(squeeze(Vv(row,1,:))), '--', ...
                    'LineWidth', 1.6, 'Color', colVar{v});
                plot(tGrid, real(squeeze(Vv(row,2,:))), '-',  ...
                    'LineWidth', 0.9, 'Color', lighten(colVar{v}));
                plot(tGrid, imag(squeeze(Vv(row,2,:))), '--', ...
                    'LineWidth', 0.9, 'Color', lighten(colVar{v}));
            end
            ylabel(rowLabel{row}, 'FontSize', 11);
            set(gca, 'XTick', 0:pi/2:Tlen, 'XTickLabel', {}); xlim([0 Tlen]);
            if row == 1
                title(sprintf(['\\nu_c^2 = %.1f, x_0 = [%g; %g]:  ', ...
                    'solid = Re, dashed = Im, thick/thin = mode 1/2'], ...
                    nuIn, x0G(1), x0G(2)), 'FontSize', 11);
                legend(hLeg, varName{1:nVar}, 'Location', 'best');
            end
        end

        % --- Tiles 3 & 4: components of w(t) = diag(e^{s t}) * V0^-1 * x0 ---
        for comp = 1:2
            nexttile; hold on; grid on;
            for v = 1:nVar
                wv = wset{v};
                lw = 1.6 - 0.7*(comp == 2);
                cc = colVar{v}; if comp == 2, cc = lighten(cc); end
                plot(tGrid, real(wv(comp,:)), '-',  'LineWidth', lw, 'Color', cc);
                plot(tGrid, imag(wv(comp,:)), '--', 'LineWidth', lw, 'Color', cc);
            end
            xlabel('Time t [rad]', 'FontSize', 11);
            ylabel(sprintf('w_%d(t) = e^{s_%d t} (V_0^{-1}x_0)_%d', comp, comp, comp), ...
                'FontSize', 11);
            set(gca, 'XTick', 0:pi/2:Tlen, ...
                'XTickLabel', strX);
            xlim([0 Tlen]);
        end

        baseNameG = sprintf('Mathieu_Factors_Case_%d_nu_%.1f_IC_%d', caseIdx, nuIn, icIdx);
        baseNameG = strrep(baseNameG, '.', 'dot');
        print(hFigG, fullfile(fDir, [baseNameG, '.svg']), '-dsvg');

%% Step H: Frequency content (FFT amplitude) of the eigenvector matrix V(t)
% The FFT runs over EXACTLY one window [0, Tlen) with the endpoint excluded,
% so the bin spacing is Omega*T/Tlen. With Tlen = 4*pi this is 0.5/rev, i.e.
% integer AND half-integer harmonics are resolved:
%   integer m      -> V is T-periodic       -> only integer harmonics
%   half-integer m -> V is anti-periodic    -> only half-integer harmonics
NfftV = 1024;
tF = linspace(0, Tlen, NfftV+1); tF(end) = [];   % endpoint excluded
dn  = T/Tlen;                                    % bin spacing in units of Omega
nAxis = (-(NfftV/2):(NfftV/2-1)) * dn;
nShow = 4.5;                                     % displayed harmonic range

% transition matrix on the FFT grid
solB = cell(1,Nz);
for k = 1:Nz
    solB{k} = ode45(@(t,x) MathieuDGL_task(t, x, D, nu02, nuC2), [t0, Tlen], I2(:,k), opts);
end
Phi_F = zeros(Nz,Nz,NfftV);
for k = 1:Nz
    Yk = deval(solB{k}, tF);
    Phi_F(:,k,:) = reshape(Yk, Nz, 1, NfftV);
end

% V(t) for each addition factor on the FFT grid (order matches Step G)
sVar  = {s_P, s_R, s_A};
VfSet = cell(1,nVar);
for v = 1:nVar
    sv = sVar{v};
    Vtmp = zeros(Nz,Nz,NfftV);
    for j = 1:NfftV
        Vtmp(:,:,j) = Phi_F(:,:,j) * V0 * ...
            diag([exp(-sv(1)*tF(j)), exp(-sv(2)*tF(j))]);
    end
    VfSet{v} = Vtmp;
end

% --- Figure: amplitude spectra, rows of V x columns (modes) ---
hFigH = figure('Name', sprintf('Case %d FFT of V(t)', caseIdx), ...
    'Color', [1 1 1], 'Position', [120 120 900 650]);
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

rowLab = {'|FFT| of V(1,:)  position row', '|FFT| of V(2,:)  velocity row'};
for row = 1:2
    for col = 1:2
        nexttile; hold on; grid on;
        hL = gobjects(1, nVar);
        for v = 1:nVar
            q = squeeze(VfSet{v}(row,col,:));
            A = abs(fftshift(fft(q)/NfftV)).';        % row vector
            sel = abs(nAxis) <= nShow;
            % small x-offset per variant so coinciding lines stay visible
            hL(v) = stem(nAxis(sel) + (v-2)*0.05, A(sel), 'filled', ...
                'Color', colVar{v}, 'MarkerSize', 3, 'LineWidth', 1.2);
        end
        xlim([-nShow nShow]); xticks(-4:4);
        if row == 2
            xlabel('Harmonic order n  [\times \Omega]', 'FontSize', 11);
        else
            set(gca, 'XTickLabel', {});
        end
        ylabel(rowLab{row}, 'FontSize', 10);
        title(sprintf('mode %d', col), 'FontSize', 10);
        if row == 1 && col == 1
            legend(hL, varName{1:nVar}, 'Location', 'best');
        end
    end
end
sgtitle(sprintf('\\nu_c^2 = %.1f:  harmonic content of the periodic eigenvector V(t)', ...
    nuIn), 'FontSize', 12, 'FontWeight', 'bold');

baseNameH = sprintf('Mathieu_FFT_V_Case_%d_nu_%.1f', caseIdx, nuIn);
baseNameH = strrep(baseNameH, '.', 'dot');
print(hFigH, fullfile(fDir, [baseNameH, '.svg']), '-dsvg');

% --- Numeric output: significant harmonics of V(1,:) ---
fprintf('  Harmonic content of V(1,:) (amplitudes > 1e-4):\n');
for v = 1:nVar
    for col = 1:2
        q = squeeze(VfSet{v}(1,col,:));
        A = abs(fftshift(fft(q)/NfftV)).';
        sel = abs(nAxis) <= 4 & A > 1e-4;
        fprintf('    %-22s mode %d: ', varName{v}, col);
        fprintf('n=%+.1f: %.4f  ', [nAxis(sel); A(sel)]);
        fprintf('\n');
    end
end

    end
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
