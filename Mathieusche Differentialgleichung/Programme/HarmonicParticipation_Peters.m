% Computes normalized harmonic modal participation vs epsilon for w = 0.7/0.3
% using mode tracking (continuation) for stability. (Figure 3 Appearance)
%
% Based on:
% David A. Peters, Sydnie M. Lieb, Loren A. Ahaus
% "Interpretation of Floquet Eigenvalues and Eigenvectors for Periodic Systems"
% JOURNAL OF THE AMERICAN HELICOPTER SOCIETY 56, 032001 (2011)
% -------------------------------------------------------------------------
clear; clc; close all;

% --- Setup for Figure Saving ---
fDir = 'figureFolder'; % Folder for figures
if ~isdir(fDir) %#ok<ISDIR>
    mkdir(fDir)
end
fDirPeters = fullfile(fDir,'figureFolderPeters'); % Subfolder specific to Peters' plots
if ~isdir(fDirPeters) %#ok<ISDIR>
    mkdir(fDirPeters)
end
% Plotting style selection
K = 'ColoredLines';
% K = 'BlackLines';
useK = strcmp(K,'BlackLines');

% -----------------------------------------------------------------------
% --- Parameters and Initialization ---
% -----------------------------------------------------------------------
Omega = 1;       % Fundamental angular frequency (Normalized: Omega = 1)
T = 2*pi / Omega;  % Period of the parametric coefficient (T = 2*pi)

% Outer loop for different unperturbed frequencies (w)
w_values = [0.3, 0.5, 0.7];
for w = w_values
    w_sq = w^2;

    % --- Basis Frequency omega0 Calculation (Peters' Convention - Eq. 10) ---
    % omega0 must satisfy 0 <= omega0 <= Omega/2 for all w values.
    basis_freq_mod = mod(w, Omega);
    if basis_freq_mod > Omega/2
        omega0 = Omega - basis_freq_mod;
    else
        omega0 = basis_freq_mod;
    end

    % Filename generation for saving the plot
    pngname = strrep(sprintf('PetersHarmonicparticipatio%s_w%2.1f',K,w),'.','dot');
    pngfile = fullfile(fDirPeters,[pngname,'.png']);

    N_FFT = 4092;    % Number of points for accurate FFT computation (should be large power of 2)
    N_eps = 400;     % Number of epsilon steps for continuation

    % Adjust eps_end based on w value for cleaner plots (matching paper's figures)
    eps_end = 5.0;
    if abs(w - 0.7) < 1e-6
        eps_end = 3.5; % Axis limit for w=0.7 plot
    end

    eps_vals = linspace(0, eps_end, N_eps); % Range of epsilon
    % Harmonics to track: m=-3..+3. These are the indices of the Fourier components.
    m_range = -3:3;
    % Color map for plotting
    colors = lines(length(m_range));
    % Storage for results
    all_participation_points = cell(N_eps, 1);

    % --- Initialization for Mode Tracking (Continuation) ---
    % The stable free vibration mode starts near eta = +/-w*Omega. We track the negative root.
    eta_initial_target = -w * Omega;   % Initial guess for the Floquet exponent (Eq. 9)
    eta_mode_prev = eta_initial_target;
    x0 = eye(2); % Initial state matrix for the Transition Matrix Phi(0) = I
    % -----------------------------------------------------------------------
    % --- Main Continuation Loop ---
    % -----------------------------------------------------------------------
    for k = 1:N_eps
        epsilon = eps_vals(k);
        % Define the State-Space matrix D(t) for the ODE (Eq. 1 & 2)
        D_func = @(t) [0, 1; -(w_sq + epsilon*sin(Omega*t)), 0];

        % Solve for the Transition Matrix Phi(T)
        sol_ode = ode45(@(t, x) reshape(D_func(t) * reshape(x, 2, 2), 4, 1), [0, T], reshape(x0, 4, 1));
        Phi_T = reshape(deval(sol_ode, T), 2, 2); % Monodromy matrix

        % Floquet Analysis: Find eigenvalues (Lambda) and eigenvectors (V)
        [V, Lambda_mat] = eig(Phi_T);
        Lambda = diag(Lambda_mat);

        % Calculate the Floquet exponents: eta = log(Lambda) / T (Eq. 9)
        eta = log(Lambda) / T;

        % --- Mode Tracking (Continuation) ---
        % Find the Floquet exponent eta_mode that is closest to the one from the previous epsilon step.
        [~, mode_idx] = min(abs(eta - eta_mode_prev));
        eta_mode = eta(mode_idx);
        v_mode = V(:, mode_idx);
        eta_mode_prev = eta_mode; % Update for the next continuation step

        % --- Compute periodic eigenvector Q(t) over one period ---
        % The Floquet solution is x(t) = exp(eta*t) * Q(t), where Q(t) is T-periodic.
        t_fft = linspace(0, T, N_FFT + 1);
        t_fft(end) = [];
        Phi_t_interp = deval(sol_ode, t_fft);
        Q_t_disp = zeros(N_FFT, 1);

        for j = 1:N_FFT
            Phi_t = reshape(Phi_t_interp(:, j), 2, 2);
            % Q(t) = Phi(t) * v_mode * exp(-eta_mode * t)
            A_t = Phi_t * v_mode * exp(-eta_mode * t_fft(j));
            Q_t_disp(j) = A_t(1); % Displacement component x(t)
        end
        % --- FFT and Harmonic Extraction (Fourier coefficients, related to Eq. 17) ---
        C = fftshift(fft(Q_t_disp) / N_FFT); % fftshift puts the zero-frequency component in the center

        freq_indices = (-N_FFT/2 : N_FFT/2 - 1);
        frequencies = freq_indices / T;

        % Extract the magnitude of the desired harmonics (m*Omega)
        harmonic_magnitudes_raw = zeros(size(m_range));
        for i = 1:length(m_range)
            m = m_range(i);
            target_freq = m * (1 / T); % Target frequency is m * Omega
            [~, idx] = min(abs(frequencies - target_freq));
            harmonic_magnitudes_raw(i) = abs(C(idx));
        end

        % --- Normalization (Modal Participation phi_m - Eq. 18) ---
        total_magnitude_sum = sum(harmonic_magnitudes_raw);
        if total_magnitude_sum > 1e-12
            phi_m = harmonic_magnitudes_raw / total_magnitude_sum; % phi_m = |C_m| / sum(|C_m|)
        else
            phi_m = zeros(size(harmonic_magnitudes_raw));
        end

        % Store results
        current_eps_points = zeros(length(m_range), 3);
        current_eps_points(:, 1) = epsilon;
        current_eps_points(:, 2) = phi_m.';
        current_eps_points(:, 3) = m_range.';
        all_participation_points{k} = current_eps_points;
    end
    % -----------------------------------------------------------------------
    % --- Plotting Section ---
    % -----------------------------------------------------------------------
    figure('Color','w','Units','pixels','Position',[200 200 900 400]);
    hold on;
    % --- Title Update: Includes ODE and w for clarity ---
    ode_str = '$\ddot{x} + (w^2 + \epsilon\sin(\Omega t)) x = 0$';
    w_str = num2str(w, '%1.1f');
    new_title = ['Harmonic participation, ', ode_str, ', $w = ', w_str, '$ ($\Omega = 2\pi$ rad/s)'];
    title(new_title, 'Interpreter', 'latex');
    xlabel('$\epsilon$', 'FontSize', 14, 'Interpreter', 'latex');
    ylabel('Modal Participation', 'FontSize', 14, 'Interpreter', 'latex');
    grid on;
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
    % Combine all data points for plotting
    all_data_matrix = vertcat(all_participation_points{:});
    % Plot each harmonic (m) as a separate curve
    jump_threshold = 0.2;
    unique_m = unique(all_data_matrix(:, 3));

    % Match line colors to the harmonic index 'm'
    m_index_map = containers.Map(unique_m, 1:length(unique_m));

    for i = 1:length(unique_m)
        m_val = unique_m(i);
        idx = (all_data_matrix(:, 3) == m_val);
        eps_for_m = all_data_matrix(idx, 1);
        phi_for_m = all_data_matrix(idx, 2);
        [eps_for_m, sort_idx] = sort(eps_for_m);
        phi_for_m = phi_for_m(sort_idx);

        % Insert NaN to break the curve where the modal participation jumps
        big_jumps = find(abs(diff(phi_for_m)) > jump_threshold);
        eps_nan = eps_for_m;
        phi_nan = phi_for_m;
        for jj = flip(big_jumps')
            eps_nan = [eps_nan(1:jj); NaN; eps_nan(jj+1:end)];
            phi_nan = [phi_nan(1:jj); NaN; phi_nan(jj+1:end)];
        end

        % --- Legend Calculation: Peters' Frequency Convention (Related to Eq. 10) ---
        if m_val >= 0
            freq_normalized = omega0/Omega + m_val;
            freq_str = sprintf('$m=%+d \\rightarrow \\omega/\\Omega \\approx %.1f$', m_val, freq_normalized);
        elseif m_val < 0
            m_abs = abs(m_val);
            freq_normalized = m_abs - omega0/Omega;
            freq_str = sprintf('$m=%+d \\rightarrow \\omega/\\Omega \\approx %.1f$', m_val, freq_normalized);
        end

        % --- Plotting Style Selection ---
        line_style = '-';
        if useK % BlackLines
            line_color = 'k';
        else % ColoredLines
            line_color = colors(m_index_map(m_val),:);
        end

        % Plot the curve for harmonic m
        plot(eps_nan, phi_nan, line_style, ...
            'Color', line_color, ...
            'LineWidth', 1.5, ...
            'DisplayName', freq_str);
    end
    axis([0 eps_end 0 1.0]);
    set(gca, 'YTick', 0:0.2:1);
    % --- Add Legend ---
    if ~useK
        legend('Location', 'northeastoutside', 'Interpreter', 'latex');
    end
    % -----------------------------------------------------------------------
    % --- Annotate labels (Adjusted position for better spacing) ---
    % -----------------------------------------------------------------------
    if abs(w - 0.3) < 1e-6
        % Labels for w = 0.3
        text(0.25, 0.8, '[+0]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
        text(0.35, 0.3, '[-1]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
        text(1.3, 0.45, '[-1/+0]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
        text(1.3, 0.15, '[-2/+1]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
        text(1.7, 0.07, '[-3/+2]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
        text(4.5, 0.18, '[-1/+1]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
        text(4.5, 0.30, '[+0]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
        text(4.0, 0.07, '[-3/+3]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
        text(4.0, 0.15, '[-2/+2]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
    elseif abs(w - 0.7) < 1e-6
        % Labels for w = 0.7
        text(0.25, 0.8, '[+0]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
        text(0.50, 0.3, '[-1]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
        text(1.3, 0.45, '[-1/+0]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
        text(1.3, 0.20, '[-2/+1]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
        text(1.7, 0.07, '[-3/+2]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
        text(3.0, 0.18, '[-1/+1]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
        text(3.0, 0.30, '[+0]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
        text(3.0, 0.07, '[-3/+3]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
        text(3.0, 0.15, '[-2/+2]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
    end
    hold off;
    % Print to Png file
    print(pngfile, '-dpng')
end