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

Omega = 1;          % Fundamental frequency
T      = 2*pi;      % Period (Omega = 1 -> T = 2*pi)

% Define all (w, epsilon) combinations
w_array       = [0.7, 0.7, 0.7, 0.7, 1.0, 1.0];
epsilon_array = [0.0, 1.0, 2.0, 3.0, 0.0, 3.5];

% Colors for each case
case_colors = {'k', 'r', 'b', 'c', 'm', 'g'};

options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

for idx = 1:length(w_array)
    w       = w_array(idx);
    epsilon = epsilon_array(idx);

    % 1. Define Mathieu ODE
    mathieu_ode = @(t, y) [y(2); -(w^2 + epsilon * sin(Omega*t)) * y(1)];

    % 2. Transition matrix Phi(T)
    [~, y1] = ode45(mathieu_ode, [0 T], [1; 0], options);
    phi_1   = y1(end, :).';
    [~, y2] = ode45(mathieu_ode, [0 T], [0; 1], options);
    phi_2   = y2(end, :).';
    Phi_T   = [phi_1, phi_2];

    % 3. Eigen-problem
    [V, D] = eig(Phi_T);
    Lambda = diag(D);

    % 4. Characteristic exponents
    eta = log(Lambda) / T;

    % Special case flags
    is_peters_case = (abs(w-0.7) < 1e-12 && abs(epsilon-3.0) < 1e-12) || ...
                     (abs(w-1.0) < 1e-12 && abs(epsilon-3.5) < 1e-12);

    % 5. Select mode index
    if (abs(w-0.7) < 1e-12 && abs(epsilon-1.0) < 1e-12)
        k = 2;
    else
        k = 1;
    end

    % 6. Standard Q(t)
    t_span = [0, 4*pi];
    x0 = V(:, k);
    if x0(1) < 0
        x0 = -x0;
    end
    [t_hist, x_hist] = ode45(mathieu_ode, t_span, x0, options);
    sigma = real(eta(k));
    Q_t   = x_hist(:, 1) .* exp(-sigma * t_hist);

    % 7. Peters-style override for special cases
    if is_peters_case
        D_func   = @(t) [0, 1; -(w^2 + epsilon*sin(Omega*t)), 0];
        D_func_i = @(t, X) reshape(D_func(t) * reshape(X,2,2), 4, 1);
        sol_Phi  = ode45(D_func_i, [0, T], reshape(eye(2),4,1));
        Phi_T_mwf = reshape(deval(sol_Phi, T), 2, 2);

        [V_mwf, D_mwf] = eig(Phi_T_mwf);
        Lambda_mwf = diag(D_mwf);
        eta_mwf    = log(Lambda_mwf) / T;

        eta_target = -w*Omega;
        [~, mode_idx_mwf] = min(abs(eta_mwf - 1i*eta_target));
        eta_mode_mwf = eta_mwf(mode_idx_mwf);
        v_mode_mwf   = V_mwf(:, mode_idx_mwf);
        v_mode_mwf   = v_mode_mwf / norm(v_mode_mwf);

        N_points_mwf = 500;
        t_vec_T = linspace(0, T, N_points_mwf);
        t_vec_T(end) = [];
        Phi_t_interp = deval(sol_Phi, t_vec_T);
        Q_t_disp_T   = zeros(length(t_vec_T),1);
        for j = 1:length(t_vec_T)
            Phi_t = reshape(Phi_t_interp(:, j), 2, 2);
            A_t   = Phi_t * v_mode_mwf * exp(-eta_mode_mwf * t_vec_T(j));
            Q_t_disp_T(j) = A_t(1);
        end

        T_plot      = 4*pi;
        num_periods = ceil(T_plot / T);
        t_hist   = [];
        Q_t_plot = [];
        Q_T_real = real(Q_t_disp_T);
        for kk = 0:num_periods-1
            t_slice = t_vec_T + kk*T;
            valid   = t_slice <= T_plot;
            if any(valid)
                t_hist   = [t_hist;   t_slice(valid)'];
                Q_t_plot = [Q_t_plot; Q_T_real(valid)];
            end
        end
        Q_t = Q_t_plot;
    end

    % 8. Plot for this (w,epsilon)
    figure('Color', 'w');
    hold on;
    plot(t_hist, real(Q_t), case_colors{idx}, 'LineWidth', 1.5);
    yline(0, 'k-');
    title(sprintf('Modal Waveform Q(t): w = %.1f, \\epsilon = %.1f', w, epsilon));
    xlabel('t (radians)');
    ylabel('Q(t)');
    xlim([0, 12.5]);
    grid on;
    legend({sprintf('Mode for w=%.1f, \\epsilon=%.1f', w, epsilon)}, 'Location', 'Best');
    hold off;

    % 9. Save this figure to PNG with a unique name
    pngname = sprintf('PetersModelWaveForm_w%1.1f_eps%1.1f', w, epsilon);
    pngname = strrep(pngname, '.', 'dot');      % avoid '.' in name
    pngfile = fullfile(fDirPeters, [pngname, '.png']);
    print(gcf, pngfile, '-dpng', '-r300');      % save current figure
end
