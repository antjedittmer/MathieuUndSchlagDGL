clear; clc; close all;

%% === PARAMETERS FROM ABBILDUNG 3.56 ===
Omega   = 1;                    % Excitation frequency
T       = 2*pi/Omega;           % Period T = 2*pi
D       = 0.15;                 % Damping factor D = 0.15
nu_vals = linspace(0, 9, 800);  % Parameter range nu_0^2 = nu_c^2
x0      = eye(2);               % Fundamental initial conditions
odeOpt  = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);

% Tie-break inside a bubble (both exponents share the same Im there)
%   'analytic'  : principal branch of sqrt(discriminant)  -> +Im branch always
%                 takes the same side of sigma = -D
%   'alternate' : tie-break flipped in every second bubble -> ladder rungs
%                 open alternately to the left and to the right
bubbleSide = 'analytic'; %'alternate'; %

%% === PRE-ALLOCATION ===
nk = numel(nu_vals);
mu_all    = zeros(nk, 2);   % col 1: +Im(s) branch, col 2: -Im(s) branch
s_R_all   = zeros(nk, 2);   % same ordering, unfolded exponents
omega_all = zeros(nk, 1);   % unfolded |Im(s)|
m_add     = zeros(nk, 1);   % addition factor m in s = log(mu)/T + i*m*Omega
logBranch = zeros(nk, 1);   % sign of the principal log branch used
isBubble  = false(nk, 1);   % real multiplier pair (tongue / bubble)

% Continuation state for the unwrapped Floquet angle theta = omega*T
sig  = +1;   % current sign of the folded angle
base = 0;    % accumulated 2*pi*m offset  -> m = base/(2*pi)
seenComplex  = false;
theta_f_prev = 0;
nBubble      = 0;    % counter of real tongues (for 'alternate')

%% === MAIN COMPUTATION LOOP ===
for k = 1:nk
    nu = nu_vals(k);

    % Mathieu ODE system matrix over period T
    ode_mat = @(t, x) [0, 1; -(nu + nu*cos(Omega*t)), -2*D] * reshape(x, 2, 2);
    [~, sol_raw] = ode45(@(t, x) reshape(ode_mat(t, x), 4, 1), [0, T], reshape(x0, 4, 1), odeOpt);
    Phi_T = reshape(sol_raw(end, :), 2, 2);

    tau  = trace(Phi_T);
    dPhi = det(Phi_T);                            % = exp(-2*D*T)
    rho  = sqrt(dPhi);

    % 1. Folded Floquet angle:  cos(theta) = tau/(2*rho)
    cth = tau / (2*rho);
    isBubble(k) = abs(cth) >= 1;                  % real pair -> theta locked at 0 or pi
    theta_f = acos(min(max(cth, -1), 1));         % folded angle, in [0, pi]

    % 2. ADDITION TERM: unwrap theta by reflection at every bubble exit.
    if k > 1 && isBubble(k-1) && ~isBubble(k) && seenComplex
        theta_prev = base + sig*theta_f_prev;
        sig  = -sig;
        base = theta_prev - sig*theta_f_prev;     % keep theta continuous
    end
    if k > 1 && ~isBubble(k-1) && isBubble(k), nBubble = nBubble + 1; end
    if ~isBubble(k), seenComplex = true; end

    theta_un = base + sig*theta_f;                % unwrapped angle
    omega    = theta_un / T;                      % Im(s_R) incl. addition term
    theta_f_prev = theta_f;

    % 3. Multipliers from the analytic branch of the discriminant.
    %    sqrt() is the principal branch, so mu_A = (tau + sqrt(disc))/2 is the
    %    continuation of the multiplier with Im >= 0 through the branch points
    %    at disc = 0. No sorting by |mu| is involved.
    sq   = sqrt(tau^2 - 4*dPhi);
    mu_A = (tau + sq)/2;
    mu_B = (tau - sq)/2;

    % Order by the imaginary part of the EXPONENT: after every reflection the
    % rung +omega is carried by the other multiplier, hence the swap with sig.
    if sig > 0
        mu_1 = mu_A; mu_2 = mu_B;
    else
        mu_1 = mu_B; mu_2 = mu_A;
    end
    if isBubble(k) && strcmp(bubbleSide, 'alternate') && mod(nBubble, 2) == 0
        [mu_1, mu_2] = deal(mu_2, mu_1);          % free tie-break
    end
    if isBubble(k) && ~seenComplex                % initial overdamped range, nu -> 0
        [mu_1, mu_2] = deal(mu_2, mu_1);
    end

    mu_all(k, :)  = [mu_1, mu_2];
    s_R_all(k, :) = [log(abs(mu_1))/T + 1i*omega, ...
                     log(abs(mu_2))/T - 1i*omega];
    omega_all(k)  = omega;
    m_add(k)      = base / (2*pi);                % integer addition factor m
    logBranch(k)  = sig;
end

% Consistency check: exp(s*T) must reproduce mu
fprintf('max |exp(s*T) - mu| = %.3e\n', max(max(abs(exp(s_R_all*T) - mu_all))));

%% === PLOTTING FIG 3.56 ===
fig356 = figure('Name', 'Abbildung 3.56', 'Color', 'w');
pos0 = get(0, 'defaultFigurePosition');
fig356.Position = [pos0(1), pos0(2)-0.15*pos0(4), pos0(3)*1.7, pos0(4)*1.1];

%% --- LEFT PLOT: FLOQUET MULTIPLIERS mu_S ---
subplot(1, 2, 1);
hold on; grid on; axis equal;

% Unit circle (|mu| = 1) and damped baseline circle (|mu| = e^(-D*T))
th = linspace(0, 2*pi, 300);
plot(cos(th), sin(th), 'k:', 'LineWidth', 1, 'DisplayName', 'Unit Circle |\mu| = 1');
r_damped = exp(-D*T);
plot(r_damped*cos(th), r_damped*sin(th), 'Color', [0.4 0.2 0.6], 'LineWidth', 1.5, ...
    'DisplayName', sprintf('Damped Radius e^{-DT} (D=%.2f)', D));

scatter(real(mu_all(:,1)), imag(mu_all(:,1)), 50, nu_vals, 'd', ...
    'DisplayName', '\mu_S (+\omega branch)');
scatter(real(mu_all(:,2)), imag(mu_all(:,2)), 25, nu_vals, 'filled', ...
    'DisplayName', '\mu_S (-\omega branch)');

xline(0, 'k--', 'Alpha', 0.3, 'HandleVisibility', 'off');
yline(0, 'k--', 'Alpha', 0.3, 'HandleVisibility', 'off');
xlim([-1.3 1.3]); ylim([-1.3 1.3]);
xlabel('\Re(\mu_S)', 'FontSize', 12);
ylabel('\Im(\mu_S)', 'FontSize', 12);
title('Floquet-Multiplikatoren \mu_S', 'FontSize', 13);
cb1 = colorbar; cb1.Label.String = '\nu_0^2 = \nu_c^2';
legend('Location', 'southoutside','NumColumns',2);
legend boxoff

%% --- RIGHT PLOT: CHARACTERISTIC EXPONENTS s_R (unfolded) ---
subplot(1, 2, 2);
hold on; grid on;

% Ladder of the addition term: Im(s) = k*Omega/2
omMax = max(abs(omega_all));
for kk = -ceil(2*omMax/Omega):ceil(2*omMax/Omega)
    yline(kk*Omega/2, '-', 'Color', 0.85*ones(1,3), 'HandleVisibility', 'off');
end

scatter(real(s_R_all(:,1)), imag(s_R_all(:,1)), 50, nu_vals, 'd', ...
    'DisplayName', 's_R (+\omega branch)');
scatter(real(s_R_all(:,2)), imag(s_R_all(:,2)), 25, nu_vals, 'filled', ...
    'DisplayName', 's_R (-\omega branch)');

xline(-D, 'r--', 'LineWidth', 1.2, 'DisplayName', sprintf('\\sigma = -D = -%.2f', D));
xline(0, 'k-', 'LineWidth', 1.0, 'DisplayName', 'Instability Boundary \sigma = 0');
xlim([-0.4 0.15]); ylim(max(3.2, 1.05*omMax)*[-1 1]);
xlabel('\Re(s_R)', 'FontSize', 12);
ylabel('\Im(s_R)', 'FontSize', 12);
title('Charakteristische Exponenten s_R', 'FontSize', 13);
cb2 = colorbar; cb2.Label.String = '\nu_0^2 = \nu_c^2';
legend('Location', 'southoutside','NumColumns',2);
legend boxoff