<<<<<<< HEAD
% Combined Script: Floquet Frequency, Harmonic Participation, and Composite Waveform
% for Mathieu's Equation.
%
% Based on:
% D.A. Peters, S.M. Lieb, L.A. Ahaus (2011)
% "Interpretation of Floquet Eigenvalues and Eigenvectors for Periodic Systems"
% Journal of the American Helicopter Society 56, 032001

clear; clc; close all;

%% ------------------------------------------------------------------------
% Setup for Figure Saving
%% ------------------------------------------------------------------------
fDir = 'figureFolder'; % Folder for figures
if ~isdir(fDir) %#ok<ISDIR>
    mkdir(fDir)
end
fDirPeters = fullfile(fDir,'figureFolderPeters'); % Subfolder specific to Peters' plots
if ~isdir(fDirPeters) %#ok<ISDIR>
    mkdir(fDirPeters)
end

% Plotting style selection
K = 'ColoredLines';        % or 'BlackLines'
useK = strcmp(K,'BlackLines');

%% ------------------------------------------------------------------------
% Parameters and Initialization
%% ------------------------------------------------------------------------
Omega = 1;           % Fundamental angular frequency (normalized)
T     = 2*pi/Omega;  % Period of the parametric coefficient (T = 2*pi)

% Outer loop for different unperturbed frequencies w
w_values =[0.3, 0.5, 0.7];
D = 0.15;
for w = w_values
    w_sq = w^2;

    % Basis frequency omega0 (Peters convention, Eq. 10)
    basis_freq_mod = mod(w, Omega);
    if basis_freq_mod > Omega/2
        omega0 = Omega - basis_freq_mod;
    else
        omega0 = basis_freq_mod;
    end

    %% ====================================================================
    % PART 1: Floquet Frequency Plot
    %% ====================================================================
    pngname = strrep(sprintf('PetersFrequency%s_w%1.1f',K,w),'.','dot');
    pngfile = fullfile(fDirPeters,[pngname,'.png']);

    % Epsilon range
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

    m_range = -4:4;  % integer multiples

    % Storage of branches by m
    results_by_branch = struct();
    for m = m_range
        if m < 0
            field_name = ['m_neg_', num2str(abs(m))];
        else
            field_name = ['m_', num2str(m)];
        end
        results_by_branch.(field_name) = [];
    end

    x0 = eye(2); % Phi(0) = I

    % Floquet exponents for each epsilon, separated into branches
    for k = 1:length(eps_vals)
        epsilon = eps_vals(k);

        % State-space matrix D(t)
        %D_func = @(t) [0, 1; -(w_sq + epsilon*sin(Omega*t)), 0];
        D_func = @(t) [0, 1; -(epsilon + epsilon*sin(Omega*t)), -2*D];


        % Solve for Phi(T)
        [~, Phi_t] = ode45(@(t, x) reshape(D_func(t)*reshape(x,2,2),4,1), ...
                            [0, T], reshape(x0,4,1));
        Phi_T = reshape(Phi_t(end,:),2,2);

        % Floquet exponents
        Lambda = eig(Phi_T);
        eta    = log(Lambda) / T;

        % Imaginary parts normalized by Omega
        normalized_omega = imag(eta) / Omega;
        basis_freq_r     = normalized_omega(1);

        % Separate into branches m
        for m = m_range
            if m == 0
                branch_freq = basis_freq_r;
            else
                branch_freq = abs(m)*Omega + sign(m)*basis_freq_r;
            end
            if m < 0
                field_name = ['m_neg_', num2str(abs(m))];
            else
                field_name = ['m_', num2str(m)];
            end
            results_by_branch.(field_name) = ...
                [results_by_branch.(field_name); epsilon, branch_freq];
        end
    end

    % Plot Part 1
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
            field_name = ['m_neg_', num2str(abs(m))];
        else
            field_name = ['m_', num2str(m)];
        end
        data = results_by_branch.(field_name);

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

    % Branch labels (as in original script)
    hold on;
    text(0.05, 0.08, '[+0]', 'FontSize',10,'Interpreter','latex');
    text(0.05, 0.55, '[-1]', 'FontSize',10,'Interpreter','latex');
    text(0.05, 1.08, '[+1]', 'FontSize',10,'Interpreter','latex');
    text(0.05, 1.55, '[-2]', 'FontSize',10,'Interpreter','latex');
    text(0.05, 2.08, '[+2]', 'FontSize',10,'Interpreter','latex');
    text(0.05, 2.55, '[-3]', 'FontSize',10,'Interpreter','latex');
    text(0.05, 3.08, '[+3]', 'FontSize',10,'Interpreter','latex');
    text(0.05, 3.55, '[-4]', 'FontSize',10,'Interpreter','latex');
    text(1.2, 0.35, '[-1/+0]', 'FontSize',10,'Interpreter','latex');
    text(1.2, 1.35, '[-2/+1]', 'FontSize',10,'Interpreter','latex');
    text(1.2, 2.35, '[-3/+2]', 'FontSize',10,'Interpreter','latex');
    text(1.2, 3.35, '[-4/+3]', 'FontSize',10,'Interpreter','latex');

    if abs(w - 0.3) < 1e-3
        text(4.0, 0.20, '[+0]',      'FontSize',10,'Interpreter','latex');
        text(4.0, 1.20, '[-1/+1]',   'FontSize',10,'Interpreter','latex');
        text(4.0, 2.20, '[-2/+2]',   'FontSize',10,'Interpreter','latex');
        text(4.0, 3.20, '[-3/+3]',   'FontSize',10,'Interpreter','latex');
        text(4.0, 4.20, '[-4/+4]',   'FontSize',10,'Interpreter','latex');
    elseif abs(w - 0.7) < 1e-3
        text(3.0, 0.20, '[+0]',      'FontSize',10,'Interpreter','latex');
        text(3.0, 1.20, '[-1/+1]',   'FontSize',10,'Interpreter','latex');
        text(3.0, 2.20, '[-2/+2]',   'FontSize',10,'Interpreter','latex');
        text(3.0, 3.20, '[-3/+3]',   'FontSize',10,'Interpreter','latex');
        text(3.0, 4.20, '[-4/+4]',   'FontSize',10,'Interpreter','latex');
    end
    print(pngfile,'-dpng');

    %% ====================================================================
    % PART 2: Harmonic Modal Participation Plot (Figure 3 Appearance)
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

    % Mode tracking: start with target eta near -w
    eta_initial_target = -w*Omega;
    eta_mode_prev = eta_initial_target;

    x0 = eye(2);

    for k = 1:N_eps
        epsilon = eps_vals_harm(k);

        D_func = @(t) [0, 1; -(w_sq + epsilon*sin(Omega*t)), 0];

        sol_ode = ode45(@(t, x) reshape(D_func(t)*reshape(x,2,2),4,1), ...
                        [0, T], reshape(x0,4,1));
        Phi_T = reshape(deval(sol_ode, T),2,2);

        [V, Lambda_mat] = eig(Phi_T);
        Lambda = diag(Lambda_mat);
        eta    = log(Lambda) / T;

        [~, mode_idx] = min(abs(eta - eta_mode_prev));
        eta_mode = eta(mode_idx);
        v_mode   = V(:,mode_idx);
        eta_mode_prev = eta_mode;

        % periodic eigenvector Q(t) over one period
        t_fft = linspace(0, T, N_FFT+1);
        t_fft(end) = [];
        Phi_t_interp = deval(sol_ode, t_fft);
        Q_t_disp = zeros(N_FFT,1);
        for j = 1:N_FFT
            Phi_t = reshape(Phi_t_interp(:,j),2,2);
            A_t   = Phi_t*v_mode*exp(-eta_mode*t_fft(j));
            Q_t_disp(j) = A_t(1);   % displacement [file:52]
        end

        % FFT and harmonic extraction
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

    % Plot Part 2: harmonic participation
    figure('Color','w','Units','pixels','Position',[200 200 900 400]);
    hold on;

    ode_str = '$\ddot{x} + (w^2 + \epsilon\sin(\Omega t)) x = 0$';
    w_str   = num2str(w,'%1.1f');
    new_title = ['Harmonic participation, ', ode_str, ...
                 ', $w = ', w_str, '$ ($\Omega = 1$ rad/s)'];
    title(new_title,'Interpreter','latex');
    xlabel('$\epsilon$','FontSize',14,'Interpreter','latex');
    ylabel('Modal Participation','FontSize',14,'Interpreter','latex');
    grid on;
    set(gca,'TickLabelInterpreter','latex','FontSize',12);

    all_data_matrix = vertcat(all_participation_points{:});
    jump_threshold = 0.2;
    unique_m = unique(all_data_matrix(:,3));
    m_index_map = containers.Map(unique_m,1:length(unique_m));

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

    if abs(w - 0.3) < 1e-6
        text(0.25,0.8,'[+0]','Interpreter','latex','FontSize',12,'FontWeight','bold');
        text(0.35,0.3,'[-1]','Interpreter','latex','FontSize',12,'FontWeight','bold');
        text(1.3,0.45,'[-1/+0]','Interpreter','latex','FontSize',12,'FontWeight','bold');
        text(1.3,0.15,'[-2/+1]','Interpreter','latex','FontSize',12,'FontWeight','bold');
        text(1.7,0.07,'[-3/+2]','Interpreter','latex','FontSize',12,'FontWeight','bold');
        text(4.5,0.18,'[-1/+1]','Interpreter','latex','FontSize',12,'FontWeight','bold');
        text(4.5,0.30,'[+0]','Interpreter','latex','FontSize',12,'FontWeight','bold');
        text(4.0,0.07,'[-3/+3]','Interpreter','latex','FontSize',12,'FontWeight','bold');
        text(4.0,0.15,'[-2/+2]','Interpreter','latex','FontSize',12,'FontWeight','bold');
    elseif abs(w - 0.7) < 1e-6
        text(0.25,0.8,'[-1]','Interpreter','latex','FontSize',12,'FontWeight','bold');
        text(0.50,0.3,'[+0]','Interpreter','latex','FontSize',12,'FontWeight','bold');
        text(1.3,0.45,'[-1/+0]','Interpreter','latex','FontSize',12,'FontWeight','bold');
        text(1.3,0.20,'[-2/+1]','Interpreter','latex','FontSize',12,'FontWeight','bold');
        text(1.7,0.07,'[-3/+2]','Interpreter','latex','FontSize',12,'FontWeight','bold');
        text(2.9,0.3,'[-1/+1]','Interpreter','latex','FontSize',12,'FontWeight','bold');
        text(2.7,0.23,'[+0]','Interpreter','latex','FontSize',12,'FontWeight','bold');
        text(3.0,0.07,'[-3/+3]','Interpreter','latex','FontSize',12,'FontWeight','bold');
        text(3.0,0.15,'[-2/+2]','Interpreter','latex','FontSize',12,'FontWeight','bold');
    end
    hold off;
    print(pngfile,'-dpng');

    %% ====================================================================
    % PART 3: Composite Modal Waveform Q_comp(t) from Harmonic Participations
    %% ====================================================================
    % Choose epsilon where you want the composite modal waveform
    eps_target = 3.0;    % adapt per w if desired
    [~, k_target] = min(abs(eps_vals_harm - eps_target));
    epsilon = eps_vals_harm(k_target);

    % Retrieve harmonic participation at this epsilon
    data_target     = all_participation_points{k_target};
    phi_m_target    = data_target(:,2);
    m_vals_target   = data_target(:,3);

    % Recompute Q(t) over one period at eps_target with same mode
    D_func = @(t) [0, 1; -(w_sq + epsilon*sin(Omega*t)), 0];
    sol_ode = ode45(@(t, x) reshape(D_func(t)*reshape(x,2,2),4,1), ...
                    [0, T], reshape(x0,4,1));
    Phi_T = reshape(deval(sol_ode, T),2,2);
    [V, Lambda_mat] = eig(Phi_T);
    Lambda = diag(Lambda_mat);
    eta    = log(Lambda) / T;
    [~, mode_idx] = min(abs(eta - eta_mode_prev));
    eta_mode = eta(mode_idx);
    v_mode   = V(:,mode_idx);

    N_comp = 4096;
    t_comp = linspace(0, T, N_comp+1);
    t_comp(end) = [];
    Phi_t_interp = deval(sol_ode, t_comp);
    Q_t_disp = zeros(N_comp,1);
    for j = 1:N_comp
        Phi_t = reshape(Phi_t_interp(:,j),2,2);
        A_t   = Phi_t*v_mode*exp(-eta_mode*t_comp(j));
        Q_t_disp(j) = A_t(1);
    end

    C = fftshift(fft(Q_t_disp)/N_comp);
    freq_indices = (-N_comp/2 : N_comp/2-1);
    frequencies  = freq_indices / T;

    C_m = zeros(size(m_vals_target));
    for i_m = 1:length(m_vals_target)
        m_val = m_vals_target(i_m);
        target_freq = m_val*(1/T);
        [~, idx_freq] = min(abs(frequencies - target_freq));
        C_m(i_m) = C(idx_freq);
    end

    C_weighted = zeros(size(C));
    for i_m = 1:length(m_vals_target)
        m_val = m_vals_target(i_m);
        target_freq = m_val*(1/T);
        [~, idx_freq] = min(abs(frequencies - target_freq));
        C_weighted(idx_freq) = phi_m_target(i_m) * exp(1i*angle(C_m(i_m)));
    end

    Q_comp = ifft(ifftshift(C_weighted*N_comp));

    % Plot composite waveform over one period
    figure('Color','w');
    plot(t_comp, real(Q_comp), 'LineWidth',1.5);
    hold on;
    yline(0,'k-');
    grid on;
    xlabel('t over one period T','Interpreter','latex');
    ylabel('$Q_{\mathrm{comp}}(t)$','Interpreter','latex');
    title(sprintf('Composite Modal Waveform $Q_{\\mathrm{comp}}(t)$, $w=%.1f$, $\\epsilon=%.2f$', ...
          w, epsilon),'Interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');
    hold off;

end
=======
% Floquet Analysis of Mathieu's Equation
% Consolidation of Frequency, Participation, and Waveform analysis.

clear; clc; close all;

%% ------------------------------------------------------------------------
% Global Parameters
%% ------------------------------------------------------------------------
Omega = 1;           
T     = 2*pi/Omega;  
w_values = [0.3, 0.5, 0.7];
N_FFT  = 4096;
m_range = -4:4;      % Combined range for frequency and participation
m_range_harm = -3:3; % Subset for modal participation plots
x0 = eye(2);         % Initial condition

% Single ODE definition: x'' + (w^2 + epsilon*sin(Omega t)) x = 0
ode_mat = @(t, epsilon, w_sq) [0, 1; -(w_sq + epsilon*sin(Omega*t)), 0];

%% ========================================================================
% Loop over w-values
%% ========================================================================
for w = w_values
    w_sq = w^2;
    
    % Setup Epsilon Range
    eps_end = 3.5; if abs(w - 0.3) < 1e-3, eps_end = 5; end
    eps_vals = linspace(0, eps_end, 400);
    
    % Pre-allocate storage
    freq_results = struct();
    for m = m_range
        f_field = matlab.lang.makeValidName(['m_' num2str(m)]);
        freq_results.(f_field) = zeros(length(eps_vals), 2);
    end
    participation_data = zeros(length(eps_vals), length(m_range_harm));
    
    % Mode tracking initialization
    eta_prev = -w * Omega;

    fprintf('Processing w = %.1f...\n', w);

    for k = 1:length(eps_vals)
        epsilon = eps_vals(k);

        % 1. Solve ODE once for this epsilon
        sol = ode45(@(t, x) reshape(ode_mat(t, epsilon, w_sq)*reshape(x,2,2),4,1), [0, T], reshape(x0,4,1));
        Phi_T = reshape(deval(sol, T), 2, 2);

        % 2. Floquet Exponents and Mode Tracking
        [V, L_mat] = eig(Phi_T);
        eta_vals = log(diag(L_mat)) / T;
        [~, idx] = min(abs(eta_vals - eta_prev)); % Track consistent branch
        eta_mode = eta_vals(idx);
        v_mode   = V(:, idx);
        eta_prev = eta_mode; 

        % 3. Harmonic Content (FFT)
        t_vec = linspace(0, T, N_FFT+1); t_vec(end) = [];
        Phi_t_all = deval(sol, t_vec);
        Q_t = zeros(N_FFT, 1);
        for j = 1:N_FFT
            Phi_curr = reshape(Phi_t_all(:,j), 2, 2);
            A_t = Phi_curr * v_mode * exp(-eta_mode * t_vec(j));
            Q_t(j) = A_t(1); % Displacement component
        end
        
        C = fftshift(fft(Q_t)/N_FFT);
        freqs = (-(N_FFT/2) : (N_FFT/2-1)) / T;

        % 4. Store Branch Frequencies
        basis_imag = imag(eta_mode) / Omega;
        for m = m_range
            f_field = matlab.lang.makeValidName(['m_' num2str(m)]);
            branch_f = (m == 0) * basis_imag + (m ~= 0) * (abs(m)*Omega + sign(m)*basis_imag);
            freq_results.(f_field)(k, :) = [epsilon, branch_f];
        end

        % 5. Store Participation
        m_mags = zeros(1, length(m_range_harm));
        for i_m = 1:length(m_range_harm)
            [~, f_idx] = min(abs(freqs - m_range_harm(i_m)/T));
            m_mags(i_m) = abs(C(f_idx));
        end
        if sum(m_mags) > 1e-12
            participation_data(k, :) = m_mags / sum(m_mags);
        end
        
        % 6. Sample Waveform (at epsilon approx 3.0)
        if k == round(length(eps_vals)*0.75) 
            figure('Name', sprintf('Waveform w=%.1f', w));
            plot(t_vec, real(Q_t), 'LineWidth', 1.5);
            title(['Waveform Q(t) at $\epsilon \approx$ ', num2str(epsilon)]);
            grid on; xlabel('Time'); ylabel('Q(t)');
        end
    end

    %% --- Plotting Results ---
    % Frequency Plot
    figure('Name', sprintf('Freq w=%.1f', w)); hold on;
    f_names = fieldnames(freq_results);
    for i = 1:length(f_names)
        data = freq_results.(f_names{i});
        plot(data(:,1), data(:,2), '.', 'MarkerSize', 4);
    end
    title(['Frequency vs \epsilon for w = ', num2str(w)]);
    xlabel('\epsilon'); ylabel('\omega / \Omega'); grid on;

    % Participation Plot
    figure('Name', sprintf('Participation w=%.1f', w));
    plot(eps_vals, participation_data, 'LineWidth', 1.5);
    title(['Harmonic Participation for w = ', num2str(w)]);
    xlabel('\epsilon'); ylabel('Participation'); grid on;
    legend(cellstr(num2str(m_range_harm', 'm=%d')));
end
>>>>>>> 16342c24dd3d857b752736907883e637769b53ac
