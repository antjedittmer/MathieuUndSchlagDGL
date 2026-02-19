% Combined Harmonic Participation for Rotor-Blade Flapping
% Covers Peters' formulation for gamma = 12.0 and gamma = 9.6
clear; clc; close all;

% --- Setup for Saving ---
fDir = 'figureFolder'; % Folder for figures
if ~isdir(fDir) %#ok<ISDIR>
    mkdir(fDir)
end
fDirPeters = fullfile(fDir,'figureFolderPeters'); % Subfolder specific to Peters' plots
if ~isdir(fDirPeters) %#ok<ISDIR>
    mkdir(fDirPeters)
end

% --- NEW: Setup for Data Saving ---
dDir = 'dataFolder'; % Folder for Excel and .mat files
if ~isdir(dDir)
    mkdir(dDir)
end

% --- Global Parameters ---
p_sq   = 1.0^2;
p      = 1.0;
Omega  = 1;
T      = 2*pi / Omega;
mu_end = 3.0;
mu_no  = 400;
m_range = -4:4;
mu_vals = linspace(0, mu_end, mu_no);

% --- Array of Gamma values to process ---
gamma_list = [12.0, 9.6];

for g_idx = 1:length(gamma_list)
    gamma = gamma_list(g_idx);
    fprintf('Processing: gamma = %.1f...\n', gamma);

    participation_matrix = zeros(length(mu_vals), length(m_range));

    % --- Main Loop over advance ratio (mu) ---
    for k = 1:length(mu_vals)
        mu = mu_vals(k);

        % State-space matrix D(t)
        D_func = @(t) [0, 1;
            -( p_sq + (gamma/8)*( (4*mu/3)*cos(t) + (mu^2)*sin(2*t) ) ), ...
            -( (gamma/8)*(1 + (4*mu/3)*sin(t)) )];

        if k == 1
            eta_mode_prev = -p * gamma;
        end

        % Solve Floquet transition matrix
        x0  = eye(2);
        rhs = @(t, x) reshape(D_func(t) * reshape(x, 2, 2), 4, 1);
        opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
        sol_ode = ode45(rhs, [0, T], reshape(x0, 4, 1), opts);

        Phi_T = reshape(deval(sol_ode, T), 2, 2);
        [V, L] = eig(Phi_T);
        L = diag(L);

        % Identify correct mode
        eta_all = log(L) / T;
        [~, mode_idx] = min(abs(eta_all - eta_mode_prev));

        eta_mode = eta_all(mode_idx);
        v_mode   = V(:, mode_idx);
        eta_mode_prev = eta_mode;

        % Periodic Modulator A(t) reconstruction
        N_fft  = 1024;
        t_fft  = linspace(0, T, N_fft+1); t_fft(end) = [];
        A_t_disp = zeros(N_fft, 1);

        for j = 1:N_fft
            Phi_t = reshape(deval(sol_ode, t_fft(j)), 2, 2);
            modulator = Phi_t * v_mode * exp(-eta_mode * t_fft(j));
            A_t_disp(j) = modulator(1);
        end

        % FFT to get harmonic strengths
        C_coeffs = fftshift(fft(A_t_disp) / N_fft);
        freq_indices = (-N_fft/2 : N_fft/2 - 1);
        freqs = freq_indices / T;

        harmonic_magnitudes_raw = zeros(size(m_range));
        for i = 1:length(m_range)
            m = m_range(i);
            target_freq = m * (1/T);
            [~, idx] = min(abs(freqs - target_freq));
            harmonic_magnitudes_raw(i) = abs(C_coeffs(idx));
        end

        % Normalize
        total_sum = sum(harmonic_magnitudes_raw);
        if total_sum > 1e-12
            participation_matrix(k, :) = harmonic_magnitudes_raw / total_sum;
        end
    end

    % --- Plotting for current gamma ---
    fig = figure('Color','w','Position', [100 100 1200 600]);
    hold on;
    color_map = lines(length(m_range));
    color_map = [color_map(1:7,:); 0 0 0; 0.5 0.5 0.5];

    jump_threshold = 0.2;
    colors = lines(length(m_range));

    for i = 1:length(m_range)
        phi_for_m = participation_matrix(:, i);
        % Handle discontinuities for clean plotting
        big_jumps = find(abs(diff(phi_for_m)) > jump_threshold);
        mu_nan  = mu_vals(:);
        phi_nan = phi_for_m(:);
        for jj = flip(big_jumps')
            mu_nan  = [mu_nan(1:jj);  NaN; mu_nan(jj+1:end)];
            phi_nan = [phi_nan(1:jj); NaN; phi_nan(jj+1:end)];
        end

        m_val = m_range(i);
        plot(mu_nan, phi_nan, '-', 'Color', color_map(i,:), 'LineWidth', 2, ...
            'DisplayName', sprintf('$m = %+d$', m_val));
    end

    % --- Apply Gamma-Specific Annotations ---
    if gamma == 12.0
        text(1.51, 0.05, '$[-4/+4]$',    'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
        text(1.22, 0.09, '$[-3/+3]$',    'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
        text(1.19, 0.15, '$[-2/+2]$',    'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
        text(1.6,  0.2, '$[-1/+1]$',     'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
        text(0.15, 0.26, '$[+0]$',       'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
        text(0.6,  0.03, '$[-4/+3]$',    'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
        text(0.5,  0.07, '$[-3/+2]$',    'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
        text(0.5,  0.20, '$[-2/+1]$',    'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
        text(0.60, 0.31, '$[-1/+0]$',    'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
        text(0.21, 0.46, '$[-1]$',       'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
    elseif gamma == 9.6
        text(1.5, 0.04, '$[-4/+4]$', 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
        text(1.4, 0.08, '$[-3/+3]$', 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
        text(1.7, 0.17, '$[-2/+2]$', 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
        text(1.5, 0.23, '$[-1/+1]$', 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
        text(0.4, 0.2, '$[+0]$', 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
        text(1.9, 0.12, '$[+0]$', 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
        text(0.6, 0.45, '$[-1]$', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
        text(0.7, 0.23, '$[+1]$', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
        text(0.3, 0.07, '$[+1]$', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
        text(0.6, 0.14, '$[-2]$', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
        text(0.7, 0.09, '$[+2]$', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
    end

    % Formatting
    title(sprintf('Harmonic Participation (Rotor Flapping: $p=1.0, \\gamma=%.1f$)', gamma), ...
        'Interpreter', 'latex', 'FontSize', 14);
    xlabel('$\mu$', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('Modal Participation', 'Interpreter', 'latex', 'FontSize', 12);
    grid on; box on;
    legend('Location', 'northeastoutside', 'Interpreter', 'latex');
    axis([0 mu_end 0 0.5]);
    set(gca, 'TickLabelInterpreter', 'latex');

    % --- File Naming and Saving ---
    gStr = strrep(num2str(gamma), '.', 'p');
    pngname = sprintf('PetersRotorParticipation_p1p0_gamma%s', gStr);

    print(fig, fullfile(fDirPeters, [pngname, '.png']), '-dpng', '-r300');

    % Save Data
    data_table = [(0:length(mu_vals)-1)', mu_vals(:), participation_matrix];
    writematrix(data_table, fullfile(dDir, [pngname, '.xlsx']));
    save(fullfile(dDir, [pngname, '.mat']), 'mu_vals', 'participation_matrix', 'm_range');

    fprintf('Done with gamma = %.1f\n\n', gamma);
end