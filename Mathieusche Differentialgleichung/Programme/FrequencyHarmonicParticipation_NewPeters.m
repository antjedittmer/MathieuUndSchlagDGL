% function FrequencyHarmonic_Peters
% Floquet Frequency Plot for Mathieu's Equation (Figure 2 Appearance)
% and Modal Harmonic Participation Plot (phi_m - Normalized Fourier Coefficients)
%
% ODE: ddot(x) + (w^2 + epsilon*sin(Omega*t)) * x = 0 (Undamped Mathieu-type)
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
Omega = 1;       % Fundamental frequency (Omega = 1 rad/s)
T = 2*pi/Omega; % Period
D = 0.15;        % Damping term, D=0.15 (Defined outside the loop)

% Outer loop for different unperturbed frequencies (w)
w_values = [0.3,0.5,0.7]; % Only w = 0.3 is processed based on user's last code change
for w = w_values
    w_sq = w^2;
    % --- Basis Frequency omega0 Calculation (Peters' Convention) ---
    % omega0 must satisfy 0 <= omega0 <= Omega/2 for all w values.
    basis_freq_mod = mod(w, Omega);
    if basis_freq_mod > Omega/2
        omega0 = Omega - basis_freq_mod;
    else
    omega0 = basis_freq_mod;
    end

    % Filename generation for saving the plot
    pngname_freq = strrep(sprintf('PetersFrequency%s_w%1.1f',K,w),'.','dot');
    pngfile_freq = fullfile(fDirPeters,[pngname_freq,'.png']);
    pngname_harm = strrep(sprintf('PetersHarmonicParticipation%s_w%1.1f',K,w),'.','dot');
    pngfile_harm = fullfile(fDirPeters,[pngname_harm,'.png']);
    
    N_FFT = 4096;    % Use power of 2 for efficient FFT
    N_eps = 400;     % Number of epsilon steps for continuation

    % Determine the range of epsilon (x-axis limit)
    if abs(w - 0.3) < 0.001
        eps_end = 5.0;
        eps_no = 1000; % High resolution for w=0.3
    elseif abs(w - 0.7) < 1e-6
        eps_end = 3.5; % Axis limit for w=0.7 plot
        eps_no = 150;
    else % Case for w=0.5 and others
        eps_end = 3.5;
        eps_no = 150;
    end
    eps_vals = linspace(0, eps_end, eps_no);
    m_range = (-3:3); % Integer multiple range for plotting branches (m*Omega)
   
    % Structure to hold frequency results, organized by branch 'm'
    results_by_branch = struct();
    for m = m_range
        if m < 0
            field_name = ['m_neg_', num2str(abs(m))];
        else
            field_name = ['m_', num2str(m)];
        end
        results_by_branch.(field_name) = []; % Initialize with empty array
    end
    
    % Initialize storage for Harmonic Participation data
    all_participation_points = cell(N_eps, 1);
    
    % -----------------------------------------------------------------------
    % --- Floquet Exponent Calculation and Branch Separation ---
    % -----------------------------------------------------------------------
    x0 = eye(2); % Initial state matrix for Phi(0) = I
    
    % Correct Initialization of eta_mode_prev: eta = +/- i*w
    eta_mode_prev = 1i * w; 

    for k = 1:length(eps_vals)
        epsilon = eps_vals(k);
        
        % D_func: ddot(x) + (w^2 + epsilon*sin(t)) * x = 0
        D_func = @(t) [0, 1; -(w_sq + epsilon*sin(Omega*t)), -2*D];
        
        % Solve for the Transition Matrix Phi(T)
        sol_ode = ode45(@(t, x) reshape(D_func(t) * reshape(x, 2, 2), 4, 1), [0, T], reshape(x0, 4, 1));
        Phi_T = reshape(deval(sol_ode, T), 2, 2); % Monodromy matrix
       
        % Calculate Floquet Exponents: eta = log(Lambda) / T
        [V, Lambda_mat] = eig(Phi_T); 
        Lambda = diag(Lambda_mat);
        eta = log(Lambda) / T;
     
        % Extract the Imaginary Part (normalized frequency)
        normalized_omega = imag(eta) / Omega;
        
        % --- Separate and Store Branches (m) (Frequency Plot Data) ---
        basis_freq_r = normalized_omega(1);
        
        for m = m_range
            if m==0
                branch_freq = basis_freq_r;
            else
                branch_freq = abs(m) *Omega + sign(m) *basis_freq_r;
            end
            if m < 0
                field_name = ['m_neg_', num2str(abs(m))];
            else
                field_name = ['m_', num2str(m)];
            end
            results_by_branch.(field_name) = [results_by_branch.(field_name); epsilon, branch_freq];
        end
        
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
        
        % --- FFT and Harmonic Extraction (Modal Participation phi_m) ---
        C = fftshift(fft(Q_t_disp) / N_FFT); % Fourier coefficients
        
        harmonic_magnitudes_raw = zeros(size(m_range));
        N_center = N_FFT/2 + 1; % Index of the zero-frequency component
        
        for i = 1:length(m_range)
            m = m_range(i);
            idx = N_center + m; 
            harmonic_magnitudes_raw(i) = abs(C(idx)); % This is Q_m^T = |C_m|
        end
        
        % --- Normalization (Modal Participation phi_m) ---
        total_magnitude_sum = sum(harmonic_magnitudes_raw);
        if total_magnitude_sum > 1e-12
            phi_m = harmonic_magnitudes_raw / total_magnitude_sum; % phi_m = |C_m| / sum(|C_m|)
        else
            phi_m = zeros(size(harmonic_magnitudes_raw));
        end
        
        % Store results at each step k
        current_eps_points = zeros(length(m_range), 3);
        current_eps_points(:, 1) = epsilon;
        current_eps_points(:, 2) = phi_m.';
        current_eps_points(:, 3) = m_range.';
        all_participation_points{k} = current_eps_points;
    end
    
    % =======================================================================
    % --- PLOTTING BLOCK 1: Frequency vs. Epsilon ---
    % =======================================================================
    figure;
    hold on;
    color_map = lines;
    
    % FIX: Use single backslashes for LaTeX in title strings
    ode_str = '$\ddot{x} + (w^2 + \epsilon\sin(\Omega t)) x = 0$';
    w_str = num2str(w, '%1.1f');
    new_title = {'Frequency vs. $\epsilon$', [ode_str, ', $w = ', w_str, '$ ($\Omega = 1$ rad/s)']};
    title(new_title, 'FontSize', 16, 'Interpreter', 'latex');
    xlabel('$\epsilon$', 'FontSize', 14, 'Interpreter', 'latex');
    ylabel('Frequency ($\omega/\Omega$)', 'FontSize', 14, 'Interpreter', 'latex');
    idx = 1;
    for m = m_range
        % Get data and plot
        if m < 0
            field_name = ['m_neg_', num2str(abs(m))];
        else
            field_name = ['m_', num2str(m)];
        end
        data = results_by_branch.(field_name);
        
        if m >= 0
            freq_normalized = omega0/Omega + m;
            freq_str = sprintf('$m=%+d \\rightarrow \\omega/\\Omega \\approx %.1f$', m, freq_normalized);
        else
            freq_normalized = abs(m) - omega0/Omega;
            freq_str = sprintf('$m=%+d \\rightarrow \\omega/\\Omega \\approx %.1f$', m, freq_normalized);
        end
        
        current_color = color_map(idx, :);
        if useK == 1
            plot(data(:, 1), data(:, 2), '.', 'Color','k', 'MarkerSize', 8, 'DisplayName', freq_str);
        else
             plot(data(:, 1), data(:, 2), '.', 'Color', current_color, 'MarkerSize', 8, 'DisplayName', freq_str);
        end
        idx = idx + 1;
    end
    grid on;
    set(gca, 'TickLabelInterpreter', 'latex');
    if ~useK
        legend('Location', 'northeastoutside', 'Interpreter', 'latex');
    end
    print(pngfile_freq, '-dpng')
    
    % =======================================================================
    % --- PLOTTING BLOCK 2: Modal Participation ($\phi_m$) vs. Epsilon ---
    % =======================================================================
    figure;
    hold on;
    % FIX: Use single backslashes for LaTeX in title strings
    new_title = {'Modal Harmonic Participation vs. $\epsilon$', ...
                 [ode_str, ', $w = ', w_str, '$ ($\Omega = 1$ rad/s)']};
    title(new_title, 'FontSize', 16, 'Interpreter', 'latex');
    xlabel('$\epsilon$', 'FontSize', 14, 'Interpreter', 'latex');
    ylabel('Modal Participation ($\phi_m$)', 'FontSize', 14, 'Interpreter', 'latex');
    grid on;
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
    
    % Combine all data points for plotting
    all_data_matrix = vertcat(all_participation_points{:});
    
    % Plot each harmonic (m) as a separate curve
    unique_m = unique(all_data_matrix(:, 3));
    m_index_map = containers.Map(unique_m, 1:length(unique_m));
    
    colors = lines(length(m_range)); % Define colors once
    
    for i = 1:length(unique_m)
        m_val = unique_m(i);
        idx = (all_data_matrix(:, 3) == m_val);
        eps_for_m = all_data_matrix(idx, 1);
        phi_for_m = all_data_matrix(idx, 2);
        
        [eps_for_m, sort_idx] = sort(eps_for_m);
        phi_for_m = phi_for_m(sort_idx);
        
        % Legend Calculation (using the standard frequency approximation for legend)
        if m_val >= 0
            freq_normalized = omega0/Omega + m_val;
            freq_str = sprintf('$m=%+d \\rightarrow \\omega/\\Omega \\approx %.1f$', m_val, freq_normalized);
        else
            freq_normalized = abs(m_val) - omega0/Omega;
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
        plot(eps_for_m, phi_for_m, line_style, ...
            'Color', line_color, ...
            'LineWidth', 1.5, ...
            'DisplayName', freq_str);
    end
    
    axis([0 eps_end 0 1.0]);
    set(gca, 'YTick', 0:0.2:1);
    
    % --- Add Legend ---
    legend('Location', 'northeastoutside', 'Interpreter', 'latex');
    
    print(pngfile_harm, '-dpng')
end