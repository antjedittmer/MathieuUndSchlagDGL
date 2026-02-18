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
dDir = fullfile('dataFolder','dataFolderPeters'); % Folder for Excel and .mat files
if ~isdir(dDir) %#ok<ISDIR>
    mkdir(dDir)
end

% --- Parameters from Peters (2011) ---
p = 1.0;            % Natural flap frequency
gamma = 9.6; %12.0;       % Lock number
N_mu = 200;
mu_vals = linspace(0, 2.5, N_mu);

m_range = -3:3;

Omega = 1;       % Fundamental angular frequency (Normalized: Omega = 1)
T = 2*pi / Omega;  % Period of the parametric coefficient (T = 2*pi)

% Storage for the specific branches shown in Figure 3
% Storage for results (Pre-allocate matrix for table creation)
participation_matrix = zeros(N_mu, length(m_range));


% Figure position for changing the position
pos0 = get(0,'defaultFigurePosition');
ode_str = 'Flap-wise';

for k = 1:N_mu
    mu = mu_vals(k);

    % Corrected periodic coefficients for Rotor Flapping
    C_t = @(t) (gamma/8) * (1 + (4/3)*mu*sin(t));
    K_t = @(t) p^2 + (gamma/8) * ((4/3)*mu*cos(t) + (mu^2)*sin(2*t));

    % State-space: x' = [0 1; -K -C]x
    A_func = @(t) [0, 1; -K_t(t), -C_t(t)];

    % 1. Compute Monodromy Matrix
    [~, X] = ode45(@(t,x) reshape(A_func(t)*reshape(x,2,2),4,1), [0 2*pi], reshape(eye(2),4,1));
    Phi_T = reshape(X(end,:), 2, 2);

    % 2. Floquet Analysis
    [V, L] = eig(Phi_T);
    % Usually we track the lead flapping mode
    [~, idx] = max(abs(diag(L)));
    lambda = L(idx,idx);
    v0 = V(:,idx);

    % Exponent (Principal value)
    eta = log(lambda)/(2*pi);

    % 3. Periodic Part Q(t) = exp(-eta*t) * Phi(t) * v0
    n_pts = 256;
    t = linspace(0, 2*pi, n_pts+1); t(end) = [];
    [~, X_path] = ode45(@(t,x) reshape(A_func(t)*reshape(x,2,2),4,1), t, reshape(eye(2),4,1));

    q_vals = zeros(n_pts, 1);
    for j = 1:n_pts
        Phi_t = reshape(X_path(j,:), 2, 2);
        q_vec = exp(-eta*t(j)) * (Phi_t * v0);
        q_vals(j) = q_vec(1); % Focus on displacement (beta)
    end

    % Perform FFT to determine normalized participation of branches m
    c = fftshift(fft(q_vals))/n_pts;

    % Create a centered frequency vector for the FFT results
    % Since C is shifted (fftshift), freq_indices must cover [-N/2, N/2-1]
    freqs = (-n_pts/2 : n_pts/2 - 1);

    % Convert indices to physical frequencies (cycles per unit time)
    % Note: frequencies are in increments of 1/T
    frequencies = freqs;

    % Pre-allocate storage for magnitudes of the selected harmonic branches
    harmonic_magnitudes_raw = zeros(size(m_range));

    for idx = 1:length(m_range)
        m = m_range(idx);

        % Define the target harmonic frequency based on Peters' n*Omega
        % Target is m * (fundamental frequency), where fundamental = 1/T
        target_freq = m;

        % Locate the FFT bin closest to the theoretical harmonic branch
        [~, idxMin] = min(abs(frequencies - target_freq));

        % Extract the magnitude |c_n| as defined in Eq. 17b
        harmonic_magnitudes_raw(idx) = abs(c(idxMin));
    end

    total_magnitude_sum = sum(harmonic_magnitudes_raw);
    if total_magnitude_sum > 1e-12
        phi_m = harmonic_magnitudes_raw / total_magnitude_sum;
    else
        phi_m = zeros(size(harmonic_magnitudes_raw));
    end

    participation_matrix(k, :) = phi_m;

end

% -----------------------------------------------------------------------
% --- Plotting Section  ---
% -----------------------------------------------------------------------
aFig = figure('Color','w','Units','pixels');
aFig.Position = [pos0(1:2), 1.4*pos0(3), pos0(4)];

hold on;

new_title = ['Harmonic participations, rotor blade flapping, $p=1.0, \gamma=' num2str(gamma) '$']; % ['Harmonic participation: ', ode_str, ', $w = ', w_str, ', \Omega = 1$\,rad/s'];
title(new_title, 'Interpreter', 'latex');
xlabel('$\mu$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Modal Participation', 'FontSize', 14, 'Interpreter', 'latex');
grid on;
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);

jump_threshold = 0.2;
line_color = 'k';
line_style = '-';
useK = 0;
colors = lines;

for i = 1:length(m_range)
    m_val = m_range(i);
    phi_for_m = participation_matrix(:, i);

    % Insert NaN to break the curve where the modal participation jumps
    % this is used for w = 0.5
    big_jumps = find(abs(diff(phi_for_m)) > jump_threshold);
    eps_nan = mu_vals;
    phi_nan = phi_for_m;
    for jj = flip(big_jumps')
        eps_nan = [eps_nan(1:jj); NaN; eps_nan(jj+1:end)];
        phi_nan = [phi_nan(1:jj); NaN; phi_nan(jj+1:end)];
    end

    % Create the legend entry for each m
    % if m_val >= 0
    %     freq_normalized = omega0/Omega + m_val;
    % elseif m_val < 0
    %     m_abs = abs(m_val);
    %     freq_normalized = m_abs - omega0/Omega;
    % end
    %freq_str = sprintf('$\\omega(m=%+d)/\\Omega \\approx %.1f$', m_val, freq_normalized);
    freq_str = sprintf('$ m=%+d $', m_val);

    % Set plotting styles and plot
    if ~useK, line_color = colors(i,:); end
    plot(eps_nan, phi_nan, line_style, 'Color', line_color, 'LineWidth', 1.5, 'DisplayName', freq_str);

end


% Add Legend
legend('Location', 'northeastoutside', 'Interpreter', 'latex');

% --- File Naming and Saving ---
gStr = strrep(num2str(gamma), '.', 'p');
pngname = sprintf('PetersHarmonicparticipationRotorFlapping_p1p0_gamma%s', gStr);
pngfile = fullfile(fDirPeters,[pngname,'.png']);
print(pngfile, '-dpng')


% % --- Plotting to match JAHS 2011 Figure 3 ---
% figure('Color','w'); hold on; grid on;
% % In the paper, curves are often solid. We use thresholds to show the 'splits'.
% % The [-1] and [0] branches start separate then merge.
% plot(mu_vals, results.m1, 'k', 'LineWidth', 1.5);
% plot(mu_vals, results.p0, 'k', 'LineWidth', 1.5);
% plot(mu_vals, results.comb21, 'k--', 'LineWidth', 1.2); % Higher order sidebands
%
% % Styling
% xlabel('Advance Ratio \mu');
% ylabel('Harmonic Participation \phi_n');
% title('Replication of Peters (2011) Fig. 10');
% axis([0 2.5 0 0.6]);