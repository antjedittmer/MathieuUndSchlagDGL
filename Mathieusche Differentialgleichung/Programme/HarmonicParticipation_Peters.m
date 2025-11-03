% modal_participation.m
% Computes normalized harmonic modal participation vs epsilon for w = 0.7/0.3
% using mode tracking (continuation) for stability.

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
w_sq = 0.7^2;    % Squared natural frequency (w^2 in the ODE)
w = sqrt(w_sq);  % Natural frequency w = 0.7 or 0.3
Omega = 1;       % Fundamental angular frequency (Omega = 1 rad/s is often used)
T = 2*pi / Omega;  % Period of the parametric coefficient (T = 2*pi)

% Define base frequency omega0 for legend/title based on Peters' convention:
% The base frequency omega0 must satisfy 0 <= omega0 <= Omega/2.
% Since w=0.3 is already in this range, it is used as omega0 for both w=0.3 and w=0.7 cases.
omega0 = 0.3; 

% Filename generation for saving the plot
pngname = strrep(sprintf('PetersHarmonicparticipatio%s_w%2.1f',K,w),'.','dot');
pngfile = fullfile(fDirPeters,[pngname,'.png']);

N_FFT = 4092;    % Number of points for accurate FFT computation (should be large power of 2)
N_eps = 400;     % Number of epsilon steps for continuation
eps_vals = linspace(0, 5.0, N_eps); % Range of epsilon (parametric excitation amplitude)
% Harmonics to track: m=-3..+3. These are the indices of the Fourier components.
m_range = -3:3;
% Color map for plotting
colors = lines(length(m_range)); 
% Storage for results: cell array to hold [epsilon, phi_m, m_val] matrices for each step
all_participation_points = cell(N_eps, 1);

% --- Initialization for Mode Tracking (Continuation) ---
% The stable free vibration mode starts near eta = +/-w*Omega. We track the negative root.
eta_initial_target = -w * Omega;   % Initial guess for the Floquet exponent (eta = mu * Omega)
eta_mode_prev = eta_initial_target; % Stores the Floquet exponent from the previous step
x0 = eye(2); % Initial state matrix for the Transition Matrix Phi(0) = I
% -----------------------------------------------------------------------
% --- Main Continuation Loop ---
% -----------------------------------------------------------------------
for k = 1:N_eps
    epsilon = eps_vals(k);

    % Define the State-Space matrix D(t) for the ODE: ddot(x) + (w^2 + epsilon*sin(Omega*t)) * x = 0
    % State vector: [x; x_dot]. D(t) is the system matrix.
    D_func = @(t) [0, 1; -(w_sq + epsilon*sin(Omega*t)), 0];
    
    % Solve the 4x4 ODE for the Transition Matrix Phi(t) = [x1, x2] (where x1, x2 are 2x1 state vectors)
    % The state is reshaped from a 2x2 matrix into a 4x1 column vector for ODE45.
    sol_ode = ode45(@(t, x) reshape(D_func(t) * reshape(x, 2, 2), 4, 1), [0, T], reshape(x0, 4, 1));
    
    % Compute the Monodromy matrix Phi(T) (the State Transition Matrix at T)
    Phi_T = reshape(deval(sol_ode, T), 2, 2); 
    
    % Floquet Analysis: Find eigenvalues (Characteristic multipliers Lambda) and eigenvectors (V)
    [V, Lambda_mat] = eig(Phi_T); 
    Lambda = diag(Lambda_mat); % Get the Char. multipliers
    
    % Calculate the Floquet exponents: eta = log(Lambda) / T.
    % The real part is the growth/decay rate, the imaginary part is the frequency.
    eta = log(Lambda) / T; 
    
    % --- CRITICAL: Mode Tracking (Continuation) ---
    % Find the Floquet exponent eta_mode that is closest to the one from the previous epsilon step.
    [~, mode_idx] = min(abs(eta - eta_mode_prev));
    eta_mode = eta(mode_idx);
    v_mode = V(:, mode_idx); % Corresponding eigenvector (v_mode is Phi(0)*v_mode for a mode)
    eta_mode_prev = eta_mode; % Update for the next continuation step
    
    % --- Compute periodic eigenvector Q(t) over one period ---
    % The Floquet solution is x(t) = exp(eta*t) * Q(t), where Q(t) is T-periodic.
    % Q(t) = Phi(t) * v_mode * exp(-eta_mode * t)
    t_fft = linspace(0, T, N_FFT + 1);
    t_fft(end) = []; % Remove the last point to avoid double counting at T=0 and T=T
    Phi_t_interp = deval(sol_ode, t_fft); % Interpolate Phi(t) at FFT points
    Q_t_disp = zeros(N_FFT, 1); % Stores the modal waveform for the displacement component (x)
    
    for j = 1:N_FFT
        Phi_t = reshape(Phi_t_interp(:, j), 2, 2);
        % Calculate the full solution x(t) = Phi(t) * v_mode * exp(eta_mode * t)
        % The periodic part Q(t) is obtained by dividing by the exponential factor.
        A_t = Phi_t * v_mode * exp(-eta_mode * t_fft(j)); 
        Q_t_disp(j) = A_t(1); % Extract the displacement component x(t)
    end

    % --- FFT and Harmonic Extraction ---
    % Compute the FFT of the periodic component Q(t)
    C = fftshift(fft(Q_t_disp) / N_FFT); % fftshift puts the zero-frequency component in the center
    
    % Determine the actual frequencies corresponding to the FFT output
    freq_indices = (-N_FFT/2 : N_FFT/2 - 1);
    frequencies = freq_indices / T; 
    
    % Extract the magnitude of the desired harmonics (m*Omega)
    harmonic_magnitudes_raw = zeros(size(m_range));
    for i = 1:length(m_range)
        m = m_range(i);
        target_freq = m * (1 / T); % Target frequency is m * Omega
        % Find the closest FFT bin to the target frequency m*Omega
        [~, idx] = min(abs(frequencies - target_freq));
        harmonic_magnitudes_raw(i) = abs(C(idx));
    end
    
    % --- Normalization (Modal Participation phi_m) ---
    % Compute the modal participation factor phi_m by normalizing the magnitudes
    % phi_m = |C_m| / sum(|C_m|)
    total_magnitude_sum = sum(harmonic_magnitudes_raw);
    if total_magnitude_sum > 1e-12
        phi_m = harmonic_magnitudes_raw / total_magnitude_sum;
    else
        phi_m = zeros(size(harmonic_magnitudes_raw));
    end
    
    % Store results: [epsilon, phi_m, m_val]
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
new_title = ['Harmonic participation, ', ode_str, ', $w = ', w_str, '$'];
title(new_title, 'Interpreter', 'latex');

xlabel('$\epsilon$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Modal Participation', 'FontSize', 14, 'Interpreter', 'latex');
grid on;
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);

% Combine all data points for plotting
all_data_matrix = vertcat(all_participation_points{:});

% Plot each harmonic (m) as a separate curve
jump_threshold = 0.2; % Threshold for inserting NaN to break discontinuous lines (mode switching)
unique_m = unique(all_data_matrix(:, 3));

for i = 1:length(unique_m)
    m_val = unique_m(i);
    idx = (all_data_matrix(:, 3) == m_val);
    eps_for_m = all_data_matrix(idx, 1);
    phi_for_m = all_data_matrix(idx, 2);
    [eps_for_m, sort_idx] = sort(eps_for_m);
    phi_for_m = phi_for_m(sort_idx);
    
    % Insert NaN to break the curve where the modal participation jumps significantly
    big_jumps = find(abs(diff(phi_for_m)) > jump_threshold);
    eps_nan = eps_for_m;
    phi_nan = phi_for_m;
    for jj = flip(big_jumps')
        eps_nan = [eps_nan(1:jj); NaN; eps_nan(jj+1:end)];
        phi_nan = [phi_nan(1:jj); NaN; phi_nan(jj+1:end)];
    end
    
    % --- Legend Calculation: Peters' Frequency Convention ---
    % The harmonic indices 'm' correspond to the frequency shift from the base frequency omega0:
    % For m >= 0: omega/Omega approx (m + omega0/Omega)
    % For m < 0: omega/Omega approx (|m| - omega0/Omega) (corresponding to the complex conjugate root)
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
        line_color = colors(i,:);
    end
    
    % Plot the curve for harmonic m
    plot(eps_nan, phi_nan, line_style, ...
         'Color', line_color, ...
         'LineWidth', 1.5, ...
         'DisplayName', freq_str); % Use the detailed frequency string in the legend
end

axis([0 5.0 0 1.0]);
set(gca, 'YTick', 0:0.2:1);

% --- Add Legend ---
if ~useK
    legend('Location', 'northeastoutside', 'Interpreter', 'latex');
end

% -----------------------------------------------------------------------
% --- Annotate labels (showing combination mode interactions) ---
% -----------------------------------------------------------------------
% These labels indicate which harmonic pair is dominating the mode shape 
% at certain instability boundaries (e.g., [-1/+0] corresponds to a mode
% shape composed of the m=-1 and m=0 harmonics).
if abs(w - 0.3) < 1e-6
    % Labels for w = 0.3
    text(0.16, 0.8, '[+0]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
    text(0.27, 0.3, '[-1]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
    text(1.2, 0.45, '[-1/+0]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
    text(1.2, 0.15, '[-2/+1]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
    text(1.6, 0.07, '[-3/+2]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
    text(4.5, 0.18, '[-1/+1]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
    text(4.5, 0.30, '[+0]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
    text(4.0, 0.07, '[-3/+3]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
    text(4.0, 0.15, '[-2/+2]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
elseif abs(w - 0.7) < 1e-6
    % Labels for w = 0.7
    text(0.19, 0.8, '[-1]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
    text(0.40, 0.3, '[+0]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
    text(1.2, 0.45, '[+0/-1]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
    text(1.2, 0.20, '[+1/-2]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
    text(1.6, 0.07, '[+2/-3]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
    text(4.5, 0.18, '[+1/-1]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
    text(4.5, 0.30, '[-1]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
    text(4.0, 0.07, '[+3/-3]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
    text(4.0, 0.15, '[+2/-2]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
end
hold off;

% Print to Png file
print(pngfile, '-dpng')