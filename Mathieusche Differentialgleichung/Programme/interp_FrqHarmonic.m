% Combined Script: Floquet Frequency, Harmonic Participation,
% and Participation-Weighted Frequency (using interp1)
%
% Based on:
% D.A. Peters, S.M. Lieb, L.A. Ahaus (2011)
% "Interpretation of Floquet Eigenvalues and Eigenvectors for Periodic Systems"
% J. Amer. Helicopter Soc. 56, 032001

clear; clc; close all;

%% ------------------------------------------------------------------------
% Setup for Figure Saving
%% ------------------------------------------------------------------------
fDir = 'figureFolder';
if ~isdir(fDir) %#ok<ISDIR>
    mkdir(fDir)
end
fDirPeters = fullfile(fDir,'figureFolderPeters');
if ~isdir(fDirPeters) %#ok<ISDIR>
    mkdir(fDirPeters)
end

K    = 'ColoredLines';          % or 'BlackLines'
useK = strcmp(K,'BlackLines');

%% ------------------------------------------------------------------------
% Global parameters
%% ------------------------------------------------------------------------
Omega = 1;           % fundamental frequency
T     = 2*pi/Omega;  % period

% Unperturbed frequencies (w) to analyze
w_values = [0.3, 0.5, 0.7];

for w = w_values
    w_sq = w^2;

    % Basis frequency omega0 (Peters convention)
    basis_freq_mod = mod(w, Omega);
    if basis_freq_mod > Omega/2
        omega0 = Omega - basis_freq_mod;
    else
        omega0 = basis_freq_mod;
    end

    %% ====================================================================
    % PART 1: Floquet Frequency Plot (Figure 2 style)
    %% ====================================================================
    pngname = strrep(sprintf('PetersFrequency%s_w%1.1f',K,w),'.','dot');
    pngfile = fullfile(fDirPeters,[pngname,'.png']);

    % epsilon range for Part 1
    if abs(w - 0.3) < 1e-3
        eps_end = 5;
        eps_no  = 1000;
    elseif abs(w - 0.7) < 1e-3
        eps_end = 3.5;
        eps_no  = 150;
    else
        eps_end = 3.5;
        eps_no  = 150;
    end
    eps_vals = linspace(0, eps_end, eps_no);

    m_range = -4:4;
    results_by_branch = struct();
    for m = m_range
        if m < 0
            fname = ['m_neg_', num2str(abs(m))];
        else
            fname = ['m_', num2str(m)];
        end
        results_by_branch.(fname) = [];
    end

    x0 = eye(2);

    % Compute Floquet exponents for each epsilon and store branches
    for k = 1:length(eps_vals)
        epsilon = eps_vals(k);
        D_func = @(t) [0, 1; -(w_sq + epsilon*sin(Omega*t)), 0];

        [~, Phi_t] = ode45(@(t,x) reshape(D_func(t)*reshape(x,2,2),4,1), ...
                           [0, T], reshape(x0,4,1));
        Phi_T = reshape(Phi_t(end,:),2,2);

        Lambda = eig(Phi_T);
        eta    = log(Lambda) / T;
        normalized_omega = imag(eta) / Omega;
        basis_freq_r     = normalized_omega(1);

        for m = m_range
            if m == 0
                branch_freq = basis_freq_r;
            else
                branch_freq = abs(m)*Omega + sign(m)*basis_freq_r;
            end
            if m < 0
                fname = ['m_neg_', num2str(abs(m))];
            else
                fname = ['m_', num2str(m)];
            end
            results_by_branch.(fname) = ...
                [results_by_branch.(fname); epsilon, branch_freq];
        end
    end

    % Plot Part 1 (optional, unchanged)
    figure;
    hold on;
    color_map = lines;
    color_map = [color_map(1:7,:); 0*ones(1,3); 0.5*ones(1,3)];

    ode_str = '$\ddot{x} + (w^2 + \epsilon\sin(\Omega t)) x = 0$';
    w_str   = num2str(w,'%1.1f');
    new_title = {'Frequency vs. $\epsilon$', ...
        [ode_str, ', $w = ', w_str, '$ ($\Omega = 1$ rad/s)']};
    title(new_title,'FontSize',16,'Interpreter','latex');
    xlabel('$\epsilon$','FontSize',14,'Interpreter','latex');
    ylabel('Frequency ($\omega/\Omega$)','FontSize',14,'Interpreter','latex');

    idx_color = 1;
    for m = m_range
        if m < 0
            fname = ['m_neg_', num2str(abs(m))];
        else
            fname = ['m_', num2str(m)];
        end
        data = results_by_branch.(fname);

        if m >= 0
            freq_normalized = omega0/Omega + m;
            freq_str = sprintf('$m=%+d \\rightarrow \\omega/\\Omega \\approx %.1f$', ...
                               m, freq_normalized);
        else
            m_abs = abs(m);
            freq_normalized = m_abs - omega0/Omega;
            freq_str = sprintf('$m=%+d \\rightarrow \\omega/\\Omega \\approx %.1f$', ...
                               m, freq_normalized);
        end

        current_color = color_map(idx_color,:);
        if useK
            plot(data(:,1), data(:,2), '.', 'Color','k', 'MarkerSize',8, ...
                 'DisplayName',freq_str);
        else
            plot(data(:,1), data(:,2), '.', 'Color',current_color, ...
                 'MarkerSize',8, 'DisplayName',freq_str);
        end
        idx_color = idx_color + 1;
    end
    grid on;
    set(gca,'TickLabelInterpreter','latex');
    if ~useK
        legend('Location','northeastoutside','Interpreter','latex');
    end
    print(pngfile,'-dpng');

    %% ====================================================================
    % PART 2: Harmonic Modal Participation (Figure 3 style)
    %% ====================================================================
    pngname = strrep(sprintf('PetersHarmonicparticipatio%s_w%2.1f',K,w),'.','dot');
    pngfile = fullfile(fDirPeters,[pngname,'.png']);

    N_FFT = 4096;
    N_eps = 400;

    eps_end_harm = 5.0;
    if abs(w - 0.7) < 1e-6
        eps_end_harm = 3.5;
    end
    eps_vals_harm = linspace(0, eps_end_harm, N_eps);

    m_range_harm = -3:3;
    colors = lines(length(m_range_harm));
    all_participation_points = cell(N_eps,1);

    eta_initial_target = -w*Omega;
    eta_mode_prev      = eta_initial_target;
    x0 = eye(2);

    for k = 1:N_eps
        epsilon = eps_vals_harm(k);
        D_func = @(t) [0, 1; -(w_sq + epsilon*sin(Omega*t)), 0];

        sol_ode = ode45(@(t,x) reshape(D_func(t)*reshape(x,2,2),4,1), ...
                        [0, T], reshape(x0,4,1));
        Phi_T = reshape(deval(sol_ode, T),2,2);

        [V, Lambda_mat] = eig(Phi_T);
        Lambda = diag(Lambda_mat);
        eta    = log(Lambda) / T;

        [~, mode_idx] = min(abs(eta - eta_mode_prev));
        eta_mode = eta(mode_idx);
        v_mode   = V(:,mode_idx);
        eta_mode_prev = eta_mode;

        t_fft = linspace(0, T, N_FFT+1);
        t_fft(end) = [];
        Phi_t_interp = deval(sol_ode, t_fft);
        Q_t_disp = zeros(N_FFT,1);
        for j = 1:N_FFT
            Phi_t = reshape(Phi_t_interp(:,j),2,2);
            A_t   = Phi_t*v_mode*exp(-eta_mode*t_fft(j));
            Q_t_disp(j) = A_t(1);
        end

        C = fftshift(fft(Q_t_disp)/N_FFT);
        freq_indices = (-N_FFT/2 : N_FFT/2-1);
        frequencies  = freq_indices / T;

        harmonic_magnitudes_raw = zeros(size(m_range_harm));
        for i_m = 1:length(m_range_harm)
            m_val = m_range_harm(i_m);
            target_freq = m_val*(1/T);
            [~, idx_freq] = min(abs(frequencies - target_freq));
            harmonic_magnitudes_raw(i_m) = abs(C(idx_freq));
        end

        total_mag = sum(harmonic_magnitudes_raw);
        if total_mag > 1e-12
            phi_m = harmonic_magnitudes_raw / total_mag;
        else
            phi_m = zeros(size(harmonic_magnitudes_raw));
        end

        current_pts = zeros(length(m_range_harm),3);
        current_pts(:,1) = epsilon;
        current_pts(:,2) = phi_m.';
        current_pts(:,3) = m_range_harm.';
        all_participation_points{k} = current_pts;
    end

    % Build all_data_matrix = [eps, phi_m, m] over all eps and m
    all_data_matrix = vertcat(all_participation_points{:});

    % Plot participation (optional, unchanged)
    figure('Color','w','Units','pixels','Position',[200 200 900 400]);
    hold on;
    w_str   = num2str(w,'%1.1f');
    new_title = ['Harmonic participation, ', ode_str, ...
                 ', $w = ', w_str, '$ ($\Omega = 1$ rad/s)'];
    title(new_title,'Interpreter','latex');
    xlabel('$\epsilon$','FontSize',14,'Interpreter','latex');
    ylabel('Modal Participation','FontSize',14,'Interpreter','latex');
    grid on;
    set(gca,'TickLabelInterpreter','latex','FontSize',12);

    jump_threshold = 0.2;
    unique_m = unique(all_data_matrix(:,3));
    m_index_map = containers.Map(unique_m, 1:length(unique_m));

    for i_m = 1:length(unique_m)
        m_val = unique_m(i_m);
        idx_m = (all_data_matrix(:,3) == m_val);
        eps_for_m = all_data_matrix(idx_m,1);
        phi_for_m = all_data_matrix(idx_m,2);
        [eps_for_m, sort_idx] = sort(eps_for_m);
        phi_for_m = phi_for_m(sort_idx);

        big_jumps = find(abs(diff(phi_for_m)) > jump_threshold);
        eps_nan = eps_for_m; phi_nan = phi_for_m;
        for jj = flip(big_jumps')
            eps_nan = [eps_nan(1:jj); NaN; eps_nan(jj+1:end)];
            phi_nan = [phi_nan(1:jj); NaN; phi_nan(jj+1:end)];
        end

        if m_val >= 0
            freq_normalized = omega0/Omega + m_val;
            freq_str = sprintf('$m=%+d \\rightarrow \\omega/\\Omega \\approx %.1f$', ...
                               m_val, freq_normalized);
        else
            m_abs = abs(m_val);
            freq_normalized = m_abs - omega0/Omega;
            freq_str = sprintf('$m=%+d \\rightarrow \\omega/\\Omega \\approx %.1f$', ...
                               m_val, freq_normalized);
        end

        line_style = '-';
        if useK
            line_color = 'k';
        else
            line_color = colors(m_index_map(m_val),:);
        end

        plot(eps_nan, phi_nan, line_style, 'Color',line_color, ...
             'LineWidth',1.5, 'DisplayName',freq_str);
    end

    axis([0 eps_end_harm 0 1.0]);
    set(gca,'YTick',0:0.2:1);
    if ~useK
        legend('Location','northeastoutside','Interpreter','latex');
    end
    hold off;
    print(pngfile,'-dpng');

    %% ====================================================================
    % PART 3: Interpolate to 150-point epsilon grid and form 150x7 matrices
    %% ====================================================================
    % Use eps grid from Part 1: take, e.g., m=0 branch
    fname0   = 'm_0';
    eps_freq = results_by_branch.(fname0)(:,1);   % 150x1 eps-grid for frequencies
    n_eps    = numel(eps_freq);
    m_modes  = -3:3;
    n_modes  = numel(m_modes);

    % 150x7 frequency matrix
    freq_mat = nan(n_eps, n_modes);
    for j = 1:n_modes
        m_val = m_modes(j);
        if m_val < 0
            fname = ['m_neg_', num2str(abs(m_val))];
        else
            fname = ['m_', num2str(m_val)];
        end
        data_branch = results_by_branch.(fname); % columns: [eps, freq]
        % If eps grids are identical, this is direct; if not, use interp1:
        % freq_mat(:,j) = interp1(data_branch(:,1), data_branch(:,2), eps_freq, 'linear','extrap');
        freq_mat(:,j) = data_branch(:,2);  % direct assignment (same eps grid)
    end

    % 150x7 participation matrix using interp1 from 400-point eps grid
    phi_mat = nan(n_eps, n_modes);
    for j = 1:n_modes
        m_val = m_modes(j);
        idx_m = (all_data_matrix(:,3) == m_val);
        eps_m = all_data_matrix(idx_m, 1);
        phi_m = all_data_matrix(idx_m, 2);
        [eps_m, sort_idx] = sort(eps_m);
        phi_m = phi_m(sort_idx);
        % interp1: map from eps_m (400) to eps_freq (150)
        phi_mat(:,j) = interp1(eps_m, phi_m, eps_freq, 'linear', 'extrap');
    end

    % Element-wise multiply and sum across modes to get weighted frequency
    summands            = freq_mat .* phi_mat;
    effective_frequency = sum(summands, 2);   % 150x1

    % Plot effective frequency vs epsilon
    figure('Color','w');
    plot(eps_freq, effective_frequency, 'k-', 'LineWidth', 1.5);
    grid on;
   xlabel('$\epsilon$','Interpreter','latex');
ylabel('Effective modal frequency (weighted)','Interpreter','latex');
title(sprintf('Weighted frequency for $w = %.1f$', w), ...
      'Interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');

end
