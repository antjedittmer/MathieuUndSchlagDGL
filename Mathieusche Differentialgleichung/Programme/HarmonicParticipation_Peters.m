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

% --- NEW: Setup for Data Saving ---
dDir = 'dataFolder'; % Folder for Excel and .mat files
if ~isdir(dDir)
    mkdir(dDir)
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
w_values =[0.3, 0.5, 0.7];
for w = w_values
    w_sq = w^2;
    % --- Basis Frequency omega0 Calculation (Peters' Convention - Eq. 10) ---
    basis_freq_mod = mod(w, Omega);
    if basis_freq_mod > Omega/2
        omega0 = Omega - basis_freq_mod;
    else
        omega0 = basis_freq_mod;
    end
    % Filename generation for saving the plot
    pngname = strrep(sprintf('PetersHarmonicparticipatio%s_w%2.1f',K,w),'.','dot');
    pngfile = fullfile(fDirPeters,[pngname,'.png']);
    
    % --- NEW: Filename generation for Data ---
    dataFileName = sprintf('Peters_HarmonicParticipation_w_%1.1f', w);
    dataFileName = strrep(dataFileName, '.', 'p'); 

    N_FFT = 4092;    
    N_eps = 400;     
    eps_end = 5.0;
    if abs(w - 0.7) < 1e-6
        eps_end = 3.5; 
    end
    eps_vals = linspace(0, eps_end, N_eps)'; % Ensure column vector
    m_range = -3:3;
    colors = lines(length(m_range));
    
    % Storage for results (Pre-allocate matrix for table creation)
    participation_matrix = zeros(N_eps, length(m_range));
    
    % --- Initialization for Mode Tracking (Continuation) ---
    eta_initial_target = -w * Omega;   
    eta_mode_prev = eta_initial_target;
    x0 = eye(2); 
    % -----------------------------------------------------------------------
    % --- Main Continuation Loop ---
    % -----------------------------------------------------------------------
    for k = 1:N_eps
        epsilon = eps_vals(k);
        D_func = @(t) [0, 1; -(w_sq + epsilon*sin(Omega*t)), 0];
        sol_ode = ode45(@(t, x) reshape(D_func(t) * reshape(x, 2, 2), 4, 1), [0, T], reshape(x0, 4, 1));
        Phi_T = reshape(deval(sol_ode, T), 2, 2); 
        [V, Lambda_mat] = eig(Phi_T);
        Lambda = diag(Lambda_mat);
        eta = log(Lambda) / T;
        [~, mode_idx] = min(abs(eta - eta_mode_prev));
        eta_mode = eta(mode_idx);
        v_mode = V(:, mode_idx);
        eta_mode_prev = eta_mode; 
        
        t_fft = linspace(0, T, N_FFT + 1);
        t_fft(end) = [];
        Phi_t_interp = deval(sol_ode, t_fft);
        Q_t_disp = zeros(N_FFT, 1);
        for j = 1:N_FFT
            Phi_t = reshape(Phi_t_interp(:, j), 2, 2);
            A_t = Phi_t * v_mode * exp(-eta_mode * t_fft(j));
            Q_t_disp(j) = A_t(1); 
        end
        C = fftshift(fft(Q_t_disp) / N_FFT); 
        freq_indices = (-N_FFT/2 : N_FFT/2 - 1);
        frequencies = freq_indices / T;
        harmonic_magnitudes_raw = zeros(size(m_range));
        for i = 1:length(m_range)
            m = m_range(i);
            target_freq = m * (1 / T); 
            [~, idx] = min(abs(frequencies - target_freq));
            harmonic_magnitudes_raw(i) = abs(C(idx));
        end
        total_magnitude_sum = sum(harmonic_magnitudes_raw);
        if total_magnitude_sum > 1e-12
            phi_m = harmonic_magnitudes_raw / total_magnitude_sum;
        else
            phi_m = zeros(size(harmonic_magnitudes_raw));
        end
        
        participation_matrix(k, :) = phi_m;
    end

    % --- NEW: Create Table and Export Data ---
    branch_names = cell(1, length(m_range));
    for i = 1:length(m_range)
        m_val = m_range(i);
        if m_val < 0
            branch_names{i} = ['m_neg_', num2str(abs(m_val))];
        else
            branch_names{i} = ['m_pos_', num2str(m_val)];
        end
    end
    
    dataTable = table(eps_vals, 'VariableNames', {'Epsilon'});
    harmonicTable = array2table(participation_matrix, 'VariableNames', branch_names);
    finalData = [dataTable, harmonicTable];
    
    writetable(finalData, fullfile(dDir, [dataFileName, '.xlsx']));
    save(fullfile(dDir, [dataFileName, '.mat']), 'finalData');

    % -----------------------------------------------------------------------
    % --- Plotting Section (Unchanged) ---
    % -----------------------------------------------------------------------
    figure('Color','w','Units','pixels','Position',[200 200 900 400]);
    hold on;
    ode_str = '$\ddot{x} + (w^2 + \epsilon\sin(\Omega t)) x = 0$';
    w_str = num2str(w, '%1.1f');
    new_title = ['Harmonic participation, ', ode_str, ', $w = ', w_str, '$ ($\Omega = 2\pi$ rad/s)'];
    title(new_title, 'Interpreter', 'latex');
    xlabel('$\epsilon$', 'FontSize', 14, 'Interpreter', 'latex');
    ylabel('Modal Participation', 'FontSize', 14, 'Interpreter', 'latex');
    grid on;
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);

    jump_threshold = 0.2;
    for i = 1:length(m_range)
        m_val = m_range(i);
        phi_for_m = participation_matrix(:, i);
        
        % Insert NaN to break the curve where the modal participation jumps
        big_jumps = find(abs(diff(phi_for_m)) > jump_threshold);
        eps_nan = eps_vals;
        phi_nan = phi_for_m;
        for jj = flip(big_jumps')
            eps_nan = [eps_nan(1:jj); NaN; eps_nan(jj+1:end)];
            phi_nan = [phi_nan(1:jj); NaN; phi_nan(jj+1:end)];
        end
        
        if m_val >= 0
            freq_normalized = omega0/Omega + m_val;
            freq_str = sprintf('$m=%+d \\rightarrow \\omega/\\Omega \\approx %.1f$', m_val, freq_normalized);
        elseif m_val < 0
            m_abs = abs(m_val);
            freq_normalized = m_abs - omega0/Omega;
            freq_str = sprintf('$m=%+d \\rightarrow \\omega/\\Omega \\approx %.1f$', m_val, freq_normalized);
        end
        
        line_style = '-';
        if useK 
            line_color = 'k';
        else 
            line_color = colors(i,:);
        end
        plot(eps_nan, phi_nan, line_style, 'Color', line_color, 'LineWidth', 1.5, 'DisplayName', freq_str);
    end
    axis([0 eps_end 0 1.0]);
    set(gca, 'YTick', 0:0.2:1);
    if ~useK
        legend('Location', 'northeastoutside', 'Interpreter', 'latex');
    end
    
    % --- Annotate labels (Unchanged) ---
    if abs(w - 0.3) < 1e-6
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
    print(pngfile, '-dpng')
end