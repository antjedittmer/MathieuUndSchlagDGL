clear; clc; close all;
%% === PARAMETERS FROM ABBILDUNG 3.56 ===
Omega = 1;                  % Excitation frequency
T = 2*pi/Omega;             % Period T = 2*pi
D = 0.15;                   % Damping factor D = 0.15
nu_vals = linspace(0, 9, 800); % Parameter range nu_0^2 = nu_c^2
x0 = eye(2);                % Fundamental initial conditions

% Pre-allocate
mu_all = zeros(length(nu_vals), 2);
s_R_all = zeros(length(nu_vals), 2);

%% === MAIN COMPUTATION LOOP ===
for k = 1:length(nu_vals)
    nu = nu_vals(k);

    % Mathieu ODE system matrix over period T
    ode_mat = @(t, x) [0, 1; -(nu + nu*cos(Omega*t)), -2*D] * reshape(x, 2, 2);
    [~, sol_raw] = ode45(@(t, x) reshape(ode_mat(t, x), 4, 1), [0, T], reshape(x0, 4, 1));
    Phi_T = reshape(sol_raw(end, :), 2, 2);

    % Floquet Multipliers mu_S
    mu_vals = eig(Phi_T);
    mu_all(k, :) = mu_vals.';

    % Characteristic Exponents s_R = ln(mu) / T  (principal branch)
    s_R_all(k, :) = log(mu_vals) / T;
end

%% === BUILD FULL LADDER OF CHARACTERISTIC EXPONENTS ===
% s_R is only unique modulo i*Omega. To reproduce the multi-branch picture
% in Abbildung 3.56 (Im axis spanning +/-3.2 while T=2*pi restricts the
% principal branch to +/-0.5), replicate each exponent across sidebands
% s_R + i*n*Omega for a range of integers n.
n_branches = -3:3;
nB = numel(n_branches);
s_R_full = zeros(length(nu_vals), 2*nB);
nu_full  = zeros(length(nu_vals), 2*nB);

col = 1;
for n = n_branches
    s_R_full(:, col)   = s_R_all(:,1) + 1i*n*Omega;
    nu_full(:, col)    = nu_vals(:);
    s_R_full(:, col+1) = s_R_all(:,2) + 1i*n*Omega;
    nu_full(:, col+1)  = nu_vals(:);
    col = col + 2;
end

%% === PLOTTING FIG 3.56 ===
fig356 = figure('Name', 'Abbildung 3.56', 'Color', 'w');
pos0 = get(0, 'defaultFigurePosition');
fig356.Position = [pos0(1), pos0(2)-0.15*pos0(4), pos0(3)*1.5, pos0(4)*1.1];

%% --- LEFT PLOT: FLOQUET MULTIPLIERS mu_S ---
subplot(1, 2, 1);
hold on; grid on; axis equal;

% Unit circle (|mu| = 1) and damped baseline circle (|mu| = e^(-D*T))
th = linspace(0, 2*pi, 300);
plot(cos(th), sin(th), 'k:', 'LineWidth', 1, 'DisplayName', 'Unit Circle |\mu| = 1');
r_damped = exp(-D*T);
plot(r_damped*cos(th), r_damped*sin(th), 'Color', [0.4 0.2 0.6], 'LineWidth', 1.5, ...
    'DisplayName', sprintf('Damped Radius e^{-DT} (D=%.2f)', D));

scatter(real(mu_all(:,1)), imag(mu_all(:,1)), 12, nu_vals, 'filled');
scatter(real(mu_all(:,2)), imag(mu_all(:,2)), 12, nu_vals, 'filled');

xline(0, 'k--', 'Alpha', 0.3, 'HandleVisibility', 'off');
yline(0, 'k--', 'Alpha', 0.3, 'HandleVisibility', 'off');
xlim([-1.3 1.3]); ylim([-1.3 1.3]);
xlabel('\Re(\mu_S)', 'FontSize', 12);
ylabel('\Im(\mu_S)', 'FontSize', 12);
title('Floquet-Multiplikatoren \mu_S', 'FontSize', 13);
cb1 = colorbar; cb1.Label.String = '\nu_0^2 = \nu_c^2';
legend('Location', 'northeast');

%% --- RIGHT PLOT: CHARACTERISTIC EXPONENTS s_R (full ladder) ---
subplot(1, 2, 2);
hold on; grid on;

scatter(real(s_R_full(:)), imag(s_R_full(:)), 12, nu_full(:), 'filled');

% Vertical backbone at Re(s) = -D
xline(-D, 'r--', 'LineWidth', 1.2, 'DisplayName', sprintf('\\sigma = -D = -%.2f', D));
xline(0, 'k-', 'LineWidth', 1.0, 'DisplayName', 'Instability Boundary \sigma = 0');
xlim([-0.4 0.15]); ylim([-3.2 3.2]);
xlabel('\Re(s_R)', 'FontSize', 12);
ylabel('\Im(s_R)', 'FontSize', 12);
title('Charakteristische Exponenten s_R', 'FontSize', 13);
cb2 = colorbar; cb2.Label.String = '\nu_0^2 = \nu_c^2';
legend('Location', 'northeast');