% Harmonic participation using SchlagDGL ODE
clear; clc; close all;

%% Paths
mathieuDir = fileparts(mfilename("fullpath"));
floquetDir = fileparts(fileparts(mathieuDir));
schlagDir  = fullfile(floquetDir, "Schlagdifferentialgleichung", "Programme");
addpath(schlagDir);

%% Global settings
SW     = 0.1;
t0     = 0.0;
plotAll = 0;
loadMat = 0;

Omega   = 1;
T       = 2*pi / Omega;
mu_min  = 0.0;
mu_max  = 3.0;
mu_no   = 400;
mu_vals = linspace(mu_min, mu_max, mu_no);
m_range = -4:4;

%% Rotor selection
AuswahlInfo = {
    1, 3, '3-Blatt-Rotor, see-saw';
    2, 3, '3-Blatt-Rotor, voll gelenkig';
    3, 4, '4-Blatt-Rotor, voll gelenkig';
    4, 5, '5-Blatt-Rotor, voll gelenkig';
    5, 3, '3-Blatt-Rotor, gelenk-/lagerlos';
    6, 4, '4-Blatt-Rotor, gelenk-/lagerlos';
    7, 1, 'Einzelblattkoordinaten im rotierenden System'
};

Auswahl = 7;   % change if needed
Blatt   = AuswahlInfo{Auswahl, 2};
rotorDescription = AuswahlInfo{Auswahl, 3};
AnzGl   = 2 * Blatt;

fprintf('Running Auswahl %d — %s\n', Auswahl, rotorDescription);

%% Load parameters
Parameter = readtable('Parameter.xlsx', 'Range', 'C4:I29');
Par       = table2array(Parameter);

ebeta = Par(8,  Auswahl);
gamma = Par(13, Auswahl);
d2    = Par(17, Auswahl);
d3    = Par(18, Auswahl);
d4    = Par(19, Auswahl);
nu0   = Par(20, Auswahl);

%% Allocation
participation_matrix = zeros(length(mu_vals), length(m_range));
eta_selected         = zeros(length(mu_vals), 1);

%% ODE settings
opts = odeset('RelTol',1e-8, 'AbsTol',1e-10, 'MaxStep',1e-3);

%% Loop over mu
for k = 1:length(mu_vals)
    mu = mu_vals(k);

    X0 = eye(AnzGl);
    sol_ode = cell(1, AnzGl);

    for j = 1:AnzGl
        sol_ode{j} = ode45( ...
            @(psi,x) SchlagDGL(psi, x, gamma, d2, d3, d4, mu, ebeta, nu0, Blatt), ...
            [0, T], X0(:,j), opts);
    end

    Phi_T = zeros(AnzGl, AnzGl);
    for j = 1:AnzGl
        Phi_T(:,j) = deval(sol_ode{j}, T);
    end

    [V, L] = eig(Phi_T);
    lambda = diag(L);
    eta_all = log(lambda) / T;

    if k == 1
        [~, mode_idx] = max(imag(eta_all));
    else
        [~, mode_idx] = min(abs(eta_all - eta_mode_prev));
    end

    eta_mode = eta_all(mode_idx);
    v_mode   = V(:, mode_idx);
    eta_mode_prev = eta_mode;
    eta_selected(k) = eta_mode;

    N_fft = 1024;
    t_fft = linspace(0, T, N_fft+1);
    t_fft(end) = [];

    A_t_disp = zeros(N_fft, 1);

    for it = 1:N_fft
        Phi_t = zeros(AnzGl, AnzGl);
        for j = 1:AnzGl
            Phi_t(:,j) = deval(sol_ode{j}, t_fft(it));
        end

        modulator = Phi_t * v_mode * exp(-eta_mode * t_fft(it));

        % observed coordinate
        A_t_disp(it) = modulator(1);
    end

    C_coeffs = fftshift(fft(A_t_disp) / N_fft);
    freq_indices = -N_fft/2 : N_fft/2 - 1;
    freqs = freq_indices / T;

    harmonic_magnitudes_raw = zeros(size(m_range));
    for i = 1:length(m_range)
        m = m_range(i);
        target_freq = m / T;
        [~, idxFreq] = min(abs(freqs - target_freq));
        harmonic_magnitudes_raw(i) = abs(C_coeffs(idxFreq));
    end

    total_sum = sum(harmonic_magnitudes_raw);
    if total_sum > 1e-12
        participation_matrix(k, :) = harmonic_magnitudes_raw / total_sum;
    end

    % fprintf('mu = %.4f done\n', mu);
end

%% Plot only
figure('Color','w', 'Position', [100 100 1200 600]);
hold on;

color_map = lines(length(m_range));
jump_threshold = 0.2;

for i = 1:length(m_range)
    phi_for_m = participation_matrix(:, i);

    big_jumps = find(abs(diff(phi_for_m)) > jump_threshold);
    mu_nan  = mu_vals(:);
    phi_nan = phi_for_m(:);

    for jj = flip(big_jumps')
        mu_nan  = [mu_nan(1:jj);  NaN; mu_nan(jj+1:end)];
        phi_nan = [phi_nan(1:jj); NaN; phi_nan(jj+1:end)];
    end

    plot(mu_nan, phi_nan, '-', ...
        'Color', color_map(i,:), ...
        'LineWidth', 2, ...
        'DisplayName', sprintf('$m = %+d$', m_range(i)));
end

title(sprintf('Harmonic Participation using SchlagDGL (%s)', rotorDescription), ...
    'Interpreter', 'none', 'FontSize', 14);
xlabel('$\mu$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Normalized harmonic participation', 'Interpreter', 'latex', 'FontSize', 12);
grid on;
box on;
legend('Location', 'northeastoutside', 'Interpreter', 'latex');
axis([mu_min mu_max 0 1]);
set(gca, 'TickLabelInterpreter', 'latex');