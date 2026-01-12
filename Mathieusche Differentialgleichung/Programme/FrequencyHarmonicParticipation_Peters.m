% Combined Script: Floquet Frequency Plot (Figure 2) and Harmonic Participation Plot (Figure 3)
% for Mathieu's Equation.
%
% Based on:
% David A. Peters, Sydnie M. Lieb, Loren A. Ahaus
% "Interpretation of Floquet Eigenvalues and Eigenvectors for Periodic Systems"
% JOURNAL OF THE AMERICAN HELICOPTER SOCIETY 56, 032001 (2011)

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
T = 2*pi/Omega; % Period of the parametric coefficient (T = 2*pi)
% Outer loop for different unperturbed frequencies (w)
w_values = [0.3, 0.5, 0.7];
for w = w_values
    w_sq = w^2;
    % --- Basis Frequency omega0 Calculation (Peters' Convention - Eq. 10) ---
    % The basis frequency, omega0, must satisfy 0 <= omega0 <= Omega/2.
    basis_freq_mod = mod(w, Omega);
    if basis_freq_mod > Omega/2
        omega0 = Omega - basis_freq_mod;
    else
        omega0 = basis_freq_mod;
    end

    % =====================================================================
    % --- PART 1: Floquet Frequency Plot (Figure 2 Appearance) ---
    % =====================================================================

    % Filename generation for saving the plot
    pngname = strrep(sprintf('PetersFrequency%s_w%1.1f',K,w),'.','dot');
    pngfile = fullfile(fDirPeters,[pngname,'.png']);
    % Determine the range of epsilon (x-axis limit)
    if abs(w - 0.3) < 0.001
        eps_end = 5;
        eps_no = 1000; % High resolution for w=0.3
    elseif abs(w - 0.7) < 0.001
        eps_end = 3.5; % Axis limit for w=0.7 plot
        eps_no = 150;
    else % Case for w=0.5 and others
        eps_end = 3.5;
        eps_no = 150;
    end
    eps_vals = linspace(0, eps_end, eps_no);
    m_range = (-4:4); % Integer multiple range for plotting branches (m*Omega)

    % Structure to hold results, organized by branch 'm'
    results_by_branch = struct();
    for m = m_range
        % Create a valid field name
        if m < 0
            field_name = ['m_neg_', num2str(abs(m))];
        else
            field_name = ['m_', num2str(m)];
        end
        results_by_branch.(field_name) = []; % Initialize with empty array
    end
    % -----------------------------------------------------------------------
    % --- Floquet Exponent Calculation and Branch Separation ---
    % -----------------------------------------------------------------------
    x0 = eye(2); % Initial state matrix for Phi(0) = I
    for k = 1:length(eps_vals)
        epsilon = eps_vals(k);
        % State-space form: d{x}/dt = [D(t)]{x}
        D_func = @(t) [0, 1; -(w_sq + epsilon*sin(Omega*t)), 0];
        % Solve for the Transition Matrix Phi(t)
        [~, Phi_t] = ode45(@(t, x) reshape(D_func(t) * reshape(x, 2, 2), 4, 1), [0, T], reshape(x0, 4, 1));
        Phi_T = reshape(Phi_t(end, :), 2, 2); % Monodromy Matrix at Phi(T)

        % Calculate Floquet Exponents: eta = log(Lambda) / T
        Lambda = eig(Phi_T);
        eta = log(Lambda) / T;

        % Extract the Imaginary Part (normalized frequency, omega/Omega = Im(eta)/Omega)
        normalized_omega = imag(eta) / Omega;
        % Choose one root (the 'basis frequency' / principal value)
        basis_freq_r = normalized_omega(1);

        % --- Separate and Store Branches (m) ---
        for m = m_range
            % This formula reconstructs the full frequency based on the integer m and basis frequency.
            % omega/Omega = abs(m) + sign(m) * (basis_freq_r/Omega)
            if m==0
                branch_freq = basis_freq_r;
            else
                branch_freq = abs(m) *Omega + sign(m) *basis_freq_r;
            end
            % Determine the valid field name
            if m < 0
                field_name = ['m_neg_', num2str(abs(m))];
            else
                field_name = ['m_', num2str(m)];
            end
            % Store [epsilon, frequency] for this specific branch
            results_by_branch.(field_name) = [results_by_branch.(field_name); epsilon, branch_freq];
        end
    end

    % -----------------------------------------------------------------------
    % --- Plotting (Part 1: Floquet Frequency) ---
    % -----------------------------------------------------------------------
    figure;
    hold on;
    color_map = lines;
    color_map = [color_map(1:7,:);0*ones(1,3);0.5*ones(1,3)];
    % Title Update
    ode_str = '$\dot{x} + (w^2 + \epsilon\sin(\Omega t)) x = 0$';
    w_str = num2str(w, '%1.1f');
    new_title = {'Frequency vs. $\epsilon$, ', [ode_str, ', $w = ', w_str, '$ ($\Omega = 1$ rad/s)']}; % Corrected Omega to 1
    title(new_title, 'FontSize', 16, 'Interpreter', 'latex');
    xlabel('$\epsilon$', 'FontSize', 14, 'Interpreter', 'latex');
    ylabel('Frequency ($\omega/\Omega$)', 'FontSize', 14, 'Interpreter', 'latex');
    idx = 1;
    for m = m_range
        % Determine field name
        if m < 0
            field_name = ['m_neg_', num2str(abs(m))];
        else
            field_name = ['m_', num2str(m)];
        end
        % --- Legend Calculation ---
        if m >= 0
            freq_normalized = omega0/Omega + m;
            freq_str = sprintf('$m=%+d \\rightarrow \\omega/\\Omega \\approx %.1f$', m, freq_normalized);
        elseif m < 0
            m_abs = abs(m);
            freq_normalized = m_abs - omega0/Omega;
            freq_str = sprintf('$m=%+d \\rightarrow \\omega/\\Omega \\approx %.1f$', m, freq_normalized);
        end
        % Get data for plotting
        data = results_by_branch.(field_name);
        current_color = color_map(idx, :);
        % Plotting
        if useK == 1
            plot(data(:, 1), data(:, 2), '.', 'Color','k', 'MarkerSize', 8, 'DisplayName', freq_str);
        else
             plot(data(:, 1), data(:, 2), '.', 'Color', current_color, 'MarkerSize', 8, 'DisplayName', freq_str);
        end
        idx = idx + 1;
    end
    % Set Axis limits
    grid on;
    set(gca, 'TickLabelInterpreter', 'latex');
    % Add Legend if colors are used
    if ~useK
        legend('Location', 'northeastoutside', 'Interpreter', 'latex');
    end

    % -----------------------------------------------------------------------
    % --- Branch Label Annotations (Part 1) ---
    % -----------------------------------------------------------------------
    % (Labels are kept as in the original script for replication purposes)
    hold on;
    % Left-side labels (ε -> 0)
    text(0.05, 0.08, '[+0]', 'FontSize', 10, 'Interpreter', 'latex');
    text(0.05, 0.55, '[-1]', 'FontSize', 10, 'Interpreter', 'latex');
    text(0.05, 1.08, '[+1]', 'FontSize', 10, 'Interpreter', 'latex');
    text(0.05, 1.55, '[-2]', 'FontSize', 10, 'Interpreter', 'latex');
    text(0.05, 2.08, '[+2]', 'FontSize', 10, 'Interpreter', 'latex');
    text(0.05, 2.55, '[-3]', 'FontSize', 10, 'Interpreter', 'latex');
    text(0.05, 3.08, '[+3]', 'FontSize', 10, 'Interpreter', 'latex');
    text(0.05, 3.55, '[-4]', 'FontSize', 10, 'Interpreter', 'latex');
    % Middle labels (Near the primary instability boundaries)
    text(1.2, 0.35, '[-1/+0]', 'FontSize', 10, 'Interpreter', 'latex');
    text(1.2, 1.35, '[-2/+1]', 'FontSize', 10, 'Interpreter', 'latex');
    text(1.2, 2.35, '[-3/+2]', 'FontSize', 10, 'Interpreter', 'latex');
    text(1.2, 3.35, '[-4/+3]', 'FontSize', 10, 'Interpreter', 'latex');
    % --- RIGHT-SIDE LABELS (Secondary Instability Boundaries) ---
    if abs(w - 0.3) < 0.001
        % Labels for w = 0.3
        text(4.0, 0.20, '[+0]', 'FontSize', 10, 'Interpreter', 'latex');
        text(4.0, 1.20, '[-1/+1]', 'FontSize', 10, 'Interpreter', 'latex');
        text(4.0, 2.20, '[-2/+2]', 'FontSize', 10, 'Interpreter', 'latex');
        text(4.0, 3.20, '[-3/+3]', 'FontSize', 10, 'Interpreter', 'latex');
        text(4.0, 4.20, '[-4/+4]', 'FontSize', 10, 'Interpreter', 'latex');
    elseif abs(w - 0.7) < 0.001
        % Labels for w = 0.7
        text(3.0, 0.20, '[+0]', 'FontSize', 10, 'Interpreter', 'latex');
        text(3.0, 1.20, '[-1/+1]', 'FontSize', 10, 'Interpreter', 'latex');
        text(3.0, 2.20, '[-2/+2]', 'FontSize', 10, 'Interpreter', 'latex');
        text(3.0, 3.20, '[-3/+3]', 'FontSize', 10, 'Interpreter', 'latex');
        text(3.0, 4.20, '[-4/+4]', 'FontSize', 10, 'Interpreter', 'latex');
    end
    % Print to Png file
    print(pngfile, '-dpng')

    % =====================================================================
    % --- PART 2: Harmonic Modal Participation Plot (Figure 3 Appearance) ---
    % =====================================================================

    % Filename generation for saving the plot
    pngname = strrep(sprintf('PetersHarmonicparticipatio%s_w%2.1f',K,w),'.','dot');
    pngfile = fullfile(fDirPeters,[pngname,'.png']);
    N_FFT = 4096;    % Number of points for accurate FFT computation (power of 2)
    N_eps = 400;     % Number of epsilon steps for continuation
    % Adjust eps_end based on w value for cleaner plots
    eps_end_harm = 5.0;
    if abs(w - 0.7) < 1e-6
        eps_end_harm = 3.5; % Axis limit for w=0.7 plot
    end
    eps_vals_harm = linspace(0, eps_end_harm, N_eps); % Range of epsilon
    % Harmonics to track: m=-3..+3.
    m_range_harm = -3:3;
    % Color map for plotting
    colors = lines(length(m_range_harm));
    % Storage for results
    all_participation_points = cell(N_eps, 1);
    % --- Initialization for Mode Tracking (Continuation) ---
    % Track the mode associated with the negative root (stable vibration mode)
    eta_initial_target = -w * Omega;
    eta_mode_prev = eta_initial_target;
    x0 = eye(2); % Initial state matrix for the Transition Matrix Phi(0) = I

    % -----------------------------------------------------------------------
    % --- Main Continuation Loop (Part 2) ---
    % -----------------------------------------------------------------------
    for k = 1:N_eps
        epsilon = eps_vals_harm(k);
        % Define the State-Space matrix D(t)
        D_func = @(t) [0, 1; -(w_sq + epsilon*sin(Omega*t)), 0];
        % Solve for the Transition Matrix Phi(T)
        sol_ode = ode45(@(t, x) reshape(D_func(t) * reshape(x, 2, 2), 4, 1), [0, T], reshape(x0, 4, 1));
        Phi_T = reshape(deval(sol_ode, T), 2, 2); % Monodromy matrix

        % Floquet Analysis
        [V, Lambda_mat] = eig(Phi_T);
        Lambda = diag(Lambda_mat);
        eta = log(Lambda) / T;

        % --- Mode Tracking (Continuation) ---
        [~, mode_idx] = min(abs(eta - eta_mode_prev));
        eta_mode = eta(mode_idx);
        v_mode = V(:, mode_idx);
        eta_mode_prev = eta_mode; % Update for the next continuation step

        % --- Compute periodic eigenvector Q(t) over one period ---
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
        % --- FFT and Harmonic Extraction ---
        C = fftshift(fft(Q_t_disp) / N_FFT); % Fourier coefficients
        freq_indices = (-N_FFT/2 : N_FFT/2 - 1);
        frequencies = freq_indices / T;
        % Extract the magnitude of the desired harmonics (m*Omega)
        harmonic_magnitudes_raw = zeros(size(m_range_harm));
        for i = 1:length(m_range_harm)
            m = m_range_harm(i);
            target_freq = m * (1 / T); % Target frequency is m * Omega
            [~, idx] = min(abs(frequencies - target_freq));
            harmonic_magnitudes_raw(i) = abs(C(idx));
        end
        % --- Normalization (Modal Participation phi_m) ---
        total_magnitude_sum = sum(harmonic_magnitudes_raw);
        if total_magnitude_sum > 1e-12
            phi_m = harmonic_magnitudes_raw / total_magnitude_sum; % phi_m = |C_m| / sum(|C_m|)
        else
            phi_m = zeros(size(harmonic_magnitudes_raw));
        end
        % Store results
        current_eps_points = zeros(length(m_range_harm), 3);
        current_eps_points(:, 1) = epsilon;
        current_eps_points(:, 2) = phi_m.';
        current_eps_points(:, 3) = m_range_harm.';
        all_participation_points{k} = current_eps_points;
    end

    % -----------------------------------------------------------------------
    % --- Plotting Section (Part 2: Harmonic Participation) ---
    % -----------------------------------------------------------------------
    figure('Color','w','Units','pixels','Position',[200 200 900 400]);
    hold on;
    % --- Title Update ---
    ode_str = '$\ddot{x} + (w^2 + \epsilon\sin(\Omega t)) x = 0$';
    w_str = num2str(w, '%1.1f');
    new_title = ['Harmonic participation, ', ode_str, ', $w = ', w_str, '$ ($\Omega = 1$ rad/s)']; % Corrected Omega to 1
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
        % --- Legend Calculation ---
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
    axis([0 eps_end_harm 0 1.0]);
    set(gca, 'YTick', 0:0.2:1);
    % --- Add Legend ---
    if ~useK
        legend('Location', 'northeastoutside', 'Interpreter', 'latex');
    end

    % -----------------------------------------------------------------------
    % --- Annotate labels (Part 2) ---
    % -----------------------------------------------------------------------
    % (Labels are kept as in the original script for replication purposes)
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
        text(0.25, 0.8, '[-1]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
        text(0.50, 0.3, '[+0]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
        text(1.3, 0.45, '[-1/+0]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
        text(1.3, 0.20, '[-2/+1]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
        text(1.7, 0.07, '[-3/+2]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
        text(2.9, 0.3, '[-1/+1]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
        text(2.7, 0.23, '[+0]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
        text(3.0, 0.07, '[-3/+3]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
        text(3.0, 0.15, '[-2/+2]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
    end
    hold off;
    % Print to Png file
    print(pngfile, '-dpng')
end