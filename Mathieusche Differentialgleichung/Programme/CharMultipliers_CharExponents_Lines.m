clear; clc; close all;

%% === GLOBAL PARAMETERS ===
Omega = 1;
T = 2*pi/Omega;
nu_vals = linspace(0, 9, 600);
m_range = -4:4;
N_FFT = 2048;
D = 0.15;
x0 = eye(2);

%% === PRE-ALLOCATION ===
multipliers_all = zeros(length(nu_vals), 2);
exponents_all = zeros(length(nu_vals), 2);

%% === MAIN COMPUTATION LOOP ===
for k = 1:length(nu_vals)
    nu = nu_vals(k);

    % Solve Monodromy Matrix over one period T
    ode_mat = @(t, x) [0, 1; -(nu + nu*cos(Omega*t)), -2*D] * reshape(x, 2, 2);
    [~, sol_raw] = ode45(@(t, x) reshape(ode_mat(t, x), 4, 1), [0, T], reshape(x0, 4, 1));
    Phi_T = reshape(sol_raw(end, :), 2, 2);

    % 1. Extract Floquet Multipliers (Eigenvalues of Monodromy Matrix)
    mu_vals = eig(Phi_T);
    multipliers_all(k, :) = mu_vals.';

    % 2. Compute Characteristic Exponents: s_R = ln(mu) / T
    s_vals = log(mu_vals) / T;
    exponents_all(k, :) = s_vals.';
end

%% === PLOTTING (REPRODUCING FIGURE 3.56) ===
pos0 = get(0, 'defaultFigurePosition');
fig = figure('Name', 'Floquet Multipliers and Exponents (Figure 3.56 Replica)', 'Color', 'w');
fig.Position = [pos0(1), pos0(2) - 0.2*pos0(4), 2*pos0(3), 1.2*pos0(4)];

tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'loose');

% --- LEFT PLOT: Floquet Multipliers \mu_s in Complex Plane ---
nexttile;
hold on; grid on; axis equal;

% Draw Unit Circle (|mu| = 1) and Inner Reference Circle (e^(-D*T))
th = linspace(0, 2*pi, 300);
plot(cos(th), sin(th), 'k--', 'LineWidth', 0.8, 'DisplayName', 'Unit Circle |\mu| = 1');
plot(exp(-D*T)*cos(th), exp(-D*T)*sin(th), 'r:', 'LineWidth', 1, 'DisplayName', sprintf('Damped Circle e^{-D T} (D=%.2f)', D));

% Scatter plot for Multipliers
scatter(real(multipliers_all(:,1)), imag(multipliers_all(:,1)), 25, nu_vals, 'filled');
scatter(real(multipliers_all(:,2)), imag(multipliers_all(:,2)), 25, nu_vals, 'filled');

% Formatting
xlabel('Re(\mu_s)', 'FontSize', 12);
ylabel('Im(\mu_s)', 'FontSize', 12);
title('Floquet Multipliers \mu_s', 'FontSize', 13);
xlim([-1.2, 1.2]); ylim([-1.2, 1.2]);
xline(0, 'k-', 'HandleVisibility', 'off'); yline(0, 'k-', 'HandleVisibility', 'off');

cb1 = colorbar;
cb1.Label.String = 'Amplification parameter \nu_c^2';

% --- RIGHT PLOT: Characteristic Exponents s_R in Complex Plane ---
nexttile;
hold on; grid on;

% Scatter plot for Positive and Negative Exponent Branches
scatter(real(exponents_all(:,1)), imag(exponents_all(:,1)), 25, nu_vals, 'filled');
scatter(real(exponents_all(:,2)), imag(exponents_all(:,2)), 25, nu_vals, 'filled');

% Highlight Stability Boundary (\sigma = 0)
xline(0, 'r--', 'LineWidth', 1.2, 'DisplayName', 'Stability Boundary \sigma = 0');
xline(-D, 'b:', 'LineWidth', 1.0, 'DisplayName', 'Mean Damping Rate \sigma = -D');

% Formatting
xlabel('Re(s_R)', 'FontSize', 12);
ylabel('Im(s_R)', 'FontSize', 12);
title('Characteristic Exponents s_R = \sigma + i\omega', 'FontSize', 13);
xlim([-0.45, 0.15]); ylim([-3.2, 3.2]);
yline(0, 'k-', 'HandleVisibility', 'off');

cb2 = colorbar;
cb2.Label.String = 'Amplification parameter \nu_c^2';