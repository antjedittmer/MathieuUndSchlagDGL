% modal_participation.m
% Computes normalized harmonic modal participation vs epsilon for w = 0.7/0.3
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
% Switch for Color/Black-and-White 
K = 'ColoredLines'; % Options: 'ColoredLines' or 'BlackLines'
useK = strcmp(K,'BlackLines'); % Flag: true if BlackLines is selected

% -----------------------------------------------------------------------
% --- Parameters and Initialization ---
% -----------------------------------------------------------------------
w_sq = 0.3^2;    % Natural frequency squared (w=0.7 to match figure 7)
w = sqrt(w_sq);  % Natural frequency w = 0.7
Omega = 1;       % Fundamental angular frequency (rad per unit time)
T = 2*pi / Omega;  % Period of the parametric coefficient (T = 2*pi)

pngname = strrep(sprintf('PetersHarmonicparticipatio%s_w%2.1f',K,w),'.','dot');
pngfile = fullfile(fDirPeters,[pngname,'.png']);

N_FFT = 4092;    % Number of points for accurate FFT 
N_eps = 400;     % Number of epsilon steps
eps_vals = linspace(0, 5.0, N_eps); % Range of epsilon (0 to 5.0)

% Harmonics to track: m=-3..+3. Figure 7 suggests tracking a mix of positive and negative indices.
m_range = -3:3;

% Color map for plotting
colors = lines(length(m_range)); % Default MATLAB colormap

% Storage for results
all_participation_points = cell(N_eps, 1);

% --- Initialization for Mode Tracking ---
% Target mode: The stable free vibration mode starts near eta = ±w.
eta_initial_target = -w * Omega;   % (Real-valued initialization, no 'i')
eta_mode_prev = eta_initial_target;

x0 = eye(2); % Initial state matrix Phi(0) = I

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
    % --- CRITICAL: Mode Tracking (Continuation) ---
    [~, mode_idx] = min(abs(eta - eta_mode_prev));
    eta_mode = eta(mode_idx);
    v_mode = V(:, mode_idx);
    eta_mode_prev = eta_mode;
    % --- Compute periodic eigenvector Q(t) over one period ---
    t_fft = linspace(0, T, N_FFT + 1);
    t_fft(end) = [];
    Phi_t_interp = deval(sol_ode, t_fft);
    Q_t_disp = zeros(N_FFT, 1);
    for j = 1:N_FFT
        Phi_t = reshape(Phi_t_interp(:, j), 2, 2);
        A_t = Phi_t * v_mode * exp(-eta_mode * t_fft(j));
        Q_t_disp(j) = A_t(1);
    end
    % --- FFT and Harmonic Extraction ---
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
    % --- Normalization (Modal Participation phi_m) ---
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
% --- Plotting Section: Updated for Color/Black and Legend ---
% -----------------------------------------------------------------------
figure('Color','w','Units','pixels','Position',[200 200 900 400]);
hold on;
% --- UPDATE: Title reflects w ---
title(['Modal Participation, $\omega = ', num2str(w, '%1.1f'), '$'], 'Interpreter', 'latex');
xlabel('$\epsilon$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Modal Participation', 'FontSize', 14, 'Interpreter', 'latex');
grid on;
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
% Combine all data points for plotting
all_data_matrix = vertcat(all_participation_points{:});
% Plot each harmonic (m) as a separate curve, breaking lines at large jumps
jump_threshold = 0.2;
unique_m = unique(all_data_matrix(:, 3));
legend_entries = {};

for i = 1:length(unique_m)
    m_val = unique_m(i);
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

    % --- COLOR/BLACK SELECTION ---
    if useK % BlackLines
        plot_style = 'k-';
        line_color = 'k';
    else % ColoredLines
        plot_style = '-';
        line_color = colors(i,:);
    end

    plot(eps_nan, phi_nan, plot_style, 'Color', line_color, 'LineWidth', 1.5, 'DisplayName', ['m=', num2str(m_val)]);
    legend_entries{end+1} = ['m=', num2str(m_val)];
end

axis([0 5.0 0 1.0]);
set(gca, 'YTick', 0:0.2:1);

% --- ADD LEGEND ONLY IF COLORED ---
if ~useK
    legend('Location', 'northeastoutside', 'Interpreter', 'latex');
end

% -----------------------------------------------------------------------
% --- Annotate labels based on frequency ---
% -----------------------------------------------------------------------
if abs(w - 0.3) < 1e-6
    % Labels for w = 0.3
    text(0.1, 0.8, '[+0]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
    text(0.15, 0.3, '[-1]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
    text(1.2, 0.45, '[-1/+0]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
    text(1.2, 0.15, '[-2/+1]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
    text(1.6, 0.07, '[-3/+2]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
    text(4.5, 0.18, '[-1/+1]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
    text(4.5, 0.30, '[+0]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
    text(4.0, 0.07, '[-3/+3]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
    text(4.0, 0.15, '[-2/+2]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
elseif abs(w - 0.7) < 1e-6
    % Labels for w = 0.7
    text(0.1, 0.8, '[-1]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
    text(0.15, 0.3, '[+0]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
    text(1.2, 0.45, '[+0/-1]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
    text(1.2, 0.15, '[+1/-2]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
    text(1.6, 0.07, '[+2/-3]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
    text(4.5, 0.18, '[+1/-1]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
    text(4.5, 0.30, '[-1]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
    text(4.0, 0.07, '[+3/-3]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
    text(4.0, 0.15, '[+2/-2]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
end
hold off;
% Print to Png file
print(pngfile, '-dpng')