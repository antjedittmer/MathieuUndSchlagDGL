% Floquet Analysis for Mathieu's Equation: multiple (w,epsilon)
clear; clc; close all;

Omega = 1;          % Fundamental frequency used in paper
T      = 2*pi;      % Period of the system (Omega = 1 -> T = 2*pi)

% Define all (w, epsilon) combinations
w_array       = 0:0.1:9; [0.7, 0.7, 0.7, 0.7, 1.0, 1.0];
epsilon_array = w_array; %[0.0, 1.0, 2.0, 3.0, 0.0, 3.5];

% Colors for each case
case_colors = {'k', 'r', 'b', 'c', 'm', 'g'};

options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

for idx = 1:length(w_array)
    w       = w_array(idx);
    epsilon = w ;epsilon_array(idx);

    % 1. Define Mathieu ODE: x'' + [w^2 + epsilon*sin(t)] x = 0
    mathieu_ode = @(t, y) [y(2); -(w^2 + epsilon * sin(Omega*t)) * y(1)];

    % 2. Compute Transition Matrix Phi(T) using 2x2 state transition
    %    (this is equivalent to integrating the 4x4 system for Phi)
    [~, y1] = ode45(mathieu_ode, [0 T], [1; 0], options);
    phi_1   = y1(end, :).';
    [~, y2] = ode45(mathieu_ode, [0 T], [0; 1], options);
    phi_2   = y2(end, :).';
    Phi_T   = [phi_1, phi_2];

    % 3. Eigen-problem (Floquet multipliers/modes)
    [V, D] = eig(Phi_T);
    Lambda = diag(D);

    % 4. Characteristic exponents
    eta = log(Lambda) / T;   % complex Floquet exponents

    % --------------------------------------------------------------------
    % CASE A: Standard "demodulated integration" method (your original)
    % --------------------------------------------------------------------
    use_standard = true;

    % For (w,eps) = (0.7,3.0) and (1.0,3.5) we will overwrite Q_t with
    % the Peters-style periodic modal waveform.
    is_peters_case = (abs(w-0.7) < 1e-12 && abs(epsilon-3.0) < 1e-12) || ...
                     (abs(w-1.0) < 1e-12 && abs(epsilon-3.5) < 1e-12);

    % 5. Select mode index for each (w,epsilon) pair
    if (abs(w-0.7) < 1e-12 && abs(epsilon-1.0) < 1e-12)
        k = 2;      % special choice matching Fig. 14b in the paper
    else
        k = 1;      % default
    end

    % 6. Standard modal waveform Q(t) using integration + demodulation
    t_span = [0, 4*pi];
    x0 = V(:, k);
    if x0(1) < 0
        x0 = -x0;
    end
    [t_hist, x_hist] = ode45(mathieu_ode, t_span, x0, options);
    sigma = real(eta(k));
    Q_t   = x_hist(:, 1) .* exp(-sigma * t_hist);   % standard Q(t) [file:52]

    % --------------------------------------------------------------------
    % CASE B: Peters-style periodic modal waveform for special cases
    %         (w,epsilon) = (0.7,3.0) and (1.0,3.5)
    % --------------------------------------------------------------------
    if is_peters_case
        % Build transition matrix Phi(t) over one period using 4x4 system
        % State vector is vec(Phi) (4x1); d/dt vec(Phi) = (I kron D(t)) vec(Phi).
        % But for 2x2, we can just propagate Phi via D(t)*Phi.
        D_func = @(t) [0, 1; -(w^2 + epsilon*sin(Omega*t)), 0];

        % ODE for Phi as 4-vector
        D_func_i = @(t, X) reshape(D_func(t) * reshape(X,2,2), 4, 1);

        % Integrate from 0 to T
        sol_Phi = ode45(D_func_i, [0, T], reshape(eye(2),4,1));
        % Recompute Phi_T from this for consistency (optional)
        Phi_T_mwf = reshape(deval(sol_Phi, T), 2, 2);  % should match Phi_T

        % Eigenanalysis again (in case of small numerical differences)
        [V_mwf, D_mwf] = eig(Phi_T_mwf);
        Lambda_mwf = diag(D_mwf);
        eta_mwf    = log(Lambda_mwf) / T;

        % Track mode whose imaginary part is closest to -w (Peters' choice)
        eta_target = -w*Omega;
        [~, mode_idx_mwf] = min(abs(eta_mwf - 1i*eta_target));
        eta_mode_mwf = eta_mwf(mode_idx_mwf);
        v_mode_mwf   = V_mwf(:, mode_idx_mwf);
        v_mode_mwf   = v_mode_mwf / norm(v_mode_mwf);

        % Compute periodic eigenvector Q(t) over one period
        N_points_mwf = 500;
        t_vec_T = linspace(0, T, N_points_mwf);
        t_vec_T(end) = [];
        Phi_t_interp = deval(sol_Phi, t_vec_T);
        Q_t_disp_T   = zeros(length(t_vec_T),1);

        for j = 1:length(t_vec_T)
            Phi_t = reshape(Phi_t_interp(:, j), 2, 2);
            A_t   = Phi_t * v_mode_mwf * exp(-eta_mode_mwf * t_vec_T(j));
            Q_t_disp_T(j) = A_t(1);     % displacement component [file:52]
        end

        % Repeat periodic waveform to fill t_span=[0,4*pi]
        T_plot      = 4*pi;
        num_periods = ceil(T_plot / T);
        t_hist      = [];
        Q_t_plot    = [];
        Q_T_real    = real(Q_t_disp_T);

        for kk = 0:num_periods-1
            t_slice = t_vec_T + kk*T;
            valid   = t_slice <= T_plot;
            if any(valid)
                t_hist   = [t_hist;   t_slice(valid)'];
                Q_t_plot = [Q_t_plot; Q_T_real(valid)];
            end
        end

        % Overwrite Q_t and t_hist with Peters-style periodic result
        Q_t   = Q_t_plot;
        % t_hist already constructed
    end

    % --------------------------------------------------------------------
    % 7. Plot (one figure per case)
    % --------------------------------------------------------------------
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
end
