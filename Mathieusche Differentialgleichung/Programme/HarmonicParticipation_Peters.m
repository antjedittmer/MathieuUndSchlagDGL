% modal_participation_w07.m
% Computes normalized harmonic modal participation vs epsilon for w = 0.3
% using mode tracking (continuation) for stability.

clear; clc; close all;

fDir = 'figureFolder'; % Ordner Abbildungen
if ~isdir(fDir) %#ok<ISDIR>
    mkdir(fDir)
end

fDirPeters = fullfile(fDir,'figureFolderPeters');
if ~isdir(fDirPeters) %#ok<ISDIR>
    mkdir(fDirPeters)
end

K = 'BlackLines';
useK = strcmp(K,'BlackLines');

% -----------------------------------------------------------------------
% --- Parameters and Initialization ---
% -----------------------------------------------------------------------
w_sq = 0.3^2;    % Natural frequency squared (w=0.7 to match figure 7)
w = sqrt(w_sq);  % Natural frequency w = 0.7
Omega = 1;       % Fundamental angular frequency (rad per unit time)
T = 2*pi / Omega;  % Period of the parametric coefficient (T = 2*pi)

pngname = strrep(sprintf('PetersHarmonicparticipatio%s_%2.1f',K,w),'.','dot');
pngfile = fullfile(fDirPeters,[pngname,'.png']);

N_FFT = 4092;    % Number of points for accurate FFT (Good value)
N_eps = 400;     % Number of epsilon steps
eps_vals = linspace(0, 5.0, N_eps); % Range of epsilon (0 to 5.0)

% Harmonics to track: m=-3..+3. Figure 7 suggests tracking a mix of positive and negative indices.
% In Floquet analysis, the *real* participation is often phi_m + phi_{-m}.
m_range = -3:3;

% Storage for results
all_participation_points = cell(N_eps, 1);

% --- Initialization for Mode Tracking ---
% Target mode: The stable free vibration mode starts at eta = +/- i*w.
% We arbitrarily track the mode starting at eta = -i*w, which typically
% corresponds to the largest 'm=0' component (phi_0 ~ 1.0).
eta_initial_target = -1 * w * Omega;
eta_mode_prev = eta_initial_target;
x0 = eye(2); % Initial state matrix Phi(0) = I

% -----------------------------------------------------------------------
% --- Main Continuation Loop ---
% -----------------------------------------------------------------------
for k = 1:N_eps
    epsilon = eps_vals(k);

    % Time-varying state matrix D(t) for x' = D(t)*x with x = [x; x_dot]
    D_func = @(t) [0, 1; -(w_sq + epsilon*sin(Omega*t)), 0]; % Corrected: Added Omega*t

    % Integrate fundamental matrix from 0 to T
    % Note: ode45 takes longer for Phi_T. A dedicated matrix exponential solver
    % or custom integrator might be faster/more stable for Phi(T) for linear systems.
    sol_ode = ode45(@(t, x) reshape(D_func(t) * reshape(x, 2, 2), 4, 1), [0, T], reshape(x0, 4, 1));
    Phi_T = reshape(deval(sol_ode, T), 2, 2);

    % Floquet multipliers and exponents
    [V, Lambda_mat] = eig(Phi_T);
    Lambda = diag(Lambda_mat);
    eta = log(Lambda) / T; % eta = sigma + i*omega

    % --- CRITICAL: Mode Tracking (Continuation) ---
    % Find the exponent (eta) closest to the one from the previous step.
    [~, mode_idx] = min(abs(eta - eta_mode_prev));
    eta_mode = eta(mode_idx);
    v_mode = V(:, mode_idx);

    % Update the previous mode for the next iteration
    eta_mode_prev = eta_mode;

    % --- Compute periodic eigenvector Q(t) over one period ---
    t_fft = linspace(0, T, N_FFT + 1); % Use N_FFT + 1 points for a period, remove last point later
    t_fft(end) = [];
    Phi_t_interp = deval(sol_ode, t_fft); % 4 x N_FFT

    Q_t_disp = zeros(N_FFT, 1); % Displacement component Q(t)_1

    for j = 1:N_FFT
        Phi_t = reshape(Phi_t_interp(:, j), 2, 2);

        % Full solution x(t) = exp(eta*t) * Q(t). We want Q(t) = exp(-eta*t) * x(t)
        % x(t) = Phi(t) * v_mode
        A_t = Phi_t * v_mode * exp(-eta_mode * t_fft(j));

        % The participation plot often uses the magnitude of the *displacement* component
        Q_t_disp(j) = A_t(1); % Use the complex displacement component
    end

    % --- FFT and Harmonic Extraction ---
    % The FFT is performed on the periodic (complex) signal Q(t).
    C = fftshift(fft(Q_t_disp) / N_FFT);
    freq_indices = (-N_FFT/2 : N_FFT/2 - 1);
    frequencies = freq_indices / T; % Frequencies in Hz (m/T)

    harmonic_magnitudes_raw = zeros(size(m_range));
    for i = 1:length(m_range)
        m = m_range(i);
        target_freq = m * (1 / T);
        [~, idx] = min(abs(frequencies - target_freq));
        harmonic_magnitudes_raw(i) = abs(C(idx)); % Use absolute value for magnitude plot
    end

    % --- Normalization (Modal Participation phi_m) ---
    % Normalization by the sum of magnitudes
    total_magnitude_sum = sum(harmonic_magnitudes_raw);
    if total_magnitude_sum > 1e-12
        phi_m = harmonic_magnitudes_raw / total_magnitude_sum;
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
% --- Plotting Section: Ensure correct mode linking ---
% -----------------------------------------------------------------------
figure('Color','w','Units','pixels','Position',[200 200 900 400]);
hold on;
title('Modal Participation, $\omega = 0.7$', 'Interpreter', 'latex');
xlabel('$\epsilon$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Modal Participation', 'FontSize', 14, 'Interpreter', 'latex');
grid on;
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);

% Combine all data points for plotting
all_data_matrix = vertcat(all_participation_points{:});

% Plot each harmonic (m) as a separate curve, breaking lines at large jumps
jump_threshold = 0.2; % You can adjust this if needed
unique_m = unique(all_data_matrix(:, 3));
for m_val = unique_m'
    idx = (all_data_matrix(:, 3) == m_val);
    eps_for_m = all_data_matrix(idx, 1);
    phi_for_m = all_data_matrix(idx, 2);
    [eps_for_m, sort_idx] = sort(eps_for_m);
    phi_for_m = phi_for_m(sort_idx);

    % Find big jumps and insert NaN to break the curve there
    big_jumps = find(abs(diff(phi_for_m)) > jump_threshold);
    eps_nan = eps_for_m;
    phi_nan = phi_for_m;
    for jj = flip(big_jumps')
        eps_nan = [eps_nan(1:jj); NaN; eps_nan(jj+1:end)];
        phi_nan = [phi_nan(1:jj); NaN; phi_nan(jj+1:end)];
    end
    plot(eps_nan, phi_nan, 'k-', 'LineWidth', 1.5);
end

axis([0 5.0 0 1.0]); % Set axis limits
set(gca, 'YTick', 0:0.2:1); % Clean up y-ticks

% Annotate labels based on Figure 7 (Visual guide only)
text(0.1, 0.8, '[+0]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
text(0.15, 0.3, '[-1]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
text(1.2, 0.45, '[-1/+0]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
text(1.2, 0.15, '[-2/+1]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
text(1.6, 0.07, '[-3/+2]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
text(4.5, 0.18, '[-1/+1]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
text(4.5, 0.30, '[+0]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
text(4.0, 0.07, '[-3/+3]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
text(4.0, 0.15, '[-2/+2]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');

hold off;

% Print to Png file
print(pngfile, '-dpng')