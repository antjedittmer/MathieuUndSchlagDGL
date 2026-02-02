% -------------------------------------------------------------------------
% Exact Harmonic Participation: Rotor Blade Flapping
% Replicating Peters' rotor flapping harmonic-participation figure
% p = 1.0, gamma = 12
% -------------------------------------------------------------------------
clear; clc; close all;

% --- Parameters (from Peters JAHS 2011) ---
p     = 1.0;          % natural flap frequency
gamma = 12.0;         % Lock number
T     = 2*pi;         % period, Omega = 1

% mu sampling
N_mu   = 250;
mu_end = 3.0;
mu_vals = linspace(0, mu_end, N_mu);

% Storage for curves
% 1: [-1], 2: [+0], 3: [-1/+0],
% 4: [-2/+1], 5: [-2/+2], 6: [-3/+2], 7: [-3/+3]
curve_data = zeros(7, N_mu);

% Target flap frequency at mu=0 (for consistent mode tracking)
target_omega0 = sqrt(1 - (gamma/16)^2);   % ≈ 0.66

% --- Loop over mu values ---
for k = 1:N_mu

    mu = mu_vals(k);

    % Periodic coefficients (JAHS Eq. (19)-(20)):
    % C(t) = gamma/8 * (1 + 4*mu/3 * sin t)
    % K(t) = p^2     + gamma/8 * (4*mu/3 * cos t + mu^2 * sin 2t)
    C_t = @(t) (gamma/8) * (1 + (4/3)*mu*sin(t));
    K_t = @(t) p^2       + (gamma/8) * ((4/3)*mu*cos(t) + mu^2*sin(2*t));

    % State-space system x' = A(t) x, x = [beta; beta_dot]
    A_func = @(t) [0, 1; -K_t(t), -C_t(t)];

    %----------------------------------------------------------------------
    % 1. Compute monodromy matrix Phi(T) over [0, 2*pi]
    %----------------------------------------------------------------------

    rhs = @(t,x) reshape(A_func(t)*reshape(x,2,2), 4, 1);
    x0  = reshape(eye(2), 4, 1);

    opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
    [~, X] = ode45(rhs, [0 2*pi], x0, opts);
    Phi_T = reshape(X(end,:), 2, 2);

    %----------------------------------------------------------------------
    % 2. Floquet eigen-analysis – choose the flap mode
    %----------------------------------------------------------------------

    [V, L]   = eig(Phi_T);
    lam_all  = diag(L);
    eta_all  = log(lam_all) / (2*pi);      % eta = sigma + i*omega
    omegas   = imag(eta_all);

    % mode index closest to target flap frequency
    [~, idx_mode] = min(abs(omegas - target_omega0));
    lambda = lam_all(idx_mode);
    v0     = V(:,idx_mode);
    eta    = log(lambda)/(2*pi);

    %----------------------------------------------------------------------
    % 3. Construct periodic part Q(t) over one period
    %   x(t) = Q(t) exp(eta t) => Q(t) = exp(-eta t) Phi(t) v0
    %----------------------------------------------------------------------

    n_fft  = 512;                          % number of time samples
    t_span = linspace(0, 2*pi, n_fft+1);
    t_span(end) = [];                      % remove duplicate point

    [~, X_path] = ode45(rhs, t_span, x0, opts);

    Q_t = zeros(n_fft, 1);
    for j = 1:n_fft
        Phi_t = reshape(X_path(j,:), 2, 2);
        q_vec = exp(-eta*t_span(j)) * (Phi_t * v0);
        Q_t(j) = q_vec(1);                 % flapping state beta
    end

    % Remove mean (to avoid artificial DC bias)
    Q_t = Q_t - mean(Q_t);

    %----------------------------------------------------------------------
    % 4. Harmonic analysis: FFT of Q(t)
    %----------------------------------------------------------------------

    % FFT, shifted so that index "mid" corresponds to harmonic m = 0
    C_fft = fftshift(fft(Q_t))/n_fft;
    mags  = abs(C_fft);
    mid   = (n_fft/2) + 1;                 % index of m=0

    % helper: magnitude of harmonic m (integer)
    get_mag = @(m) mags(mid + m);

    total_mag = sum(mags);
    if total_mag < eps
        % degenerate case (should not occur for normal parameters)
        curve_data(:,k) = 0;
        continue;
    end

    % Individual [-1] and [0]
    m_neg1 = get_mag(-1);
    m_0    = get_mag( 0);

    % Combined branches
    comb_m1_0  = m_neg1 + m_0;           % [-1/+0]
    comb_m2_p1 = get_mag(-2) + get_mag(1);  % [-2/+1]
    comb_m2_p2 = get_mag(-2) + get_mag(2);  % [-2/+2]
    comb_m3_p2 = get_mag(-3) + get_mag(2);  % [-3/+2]
    comb_m3_p3 = get_mag(-3) + get_mag(3);  % [-3/+3]

    % Normalize by sum of all |c_n| to get participation (φ_n)[file:2]
    curve_data(1,k) = m_neg1 / total_mag;         % [-1]
    curve_data(2,k) = m_0    / total_mag;         % [+0]
    curve_data(3,k) = comb_m1_0  / total_mag;     % [-1/+0]
    curve_data(4,k) = comb_m2_p1 / total_mag;     % [-2/+1]
    curve_data(5,k) = comb_m2_p2 / total_mag;     % [-2/+2]
    curve_data(6,k) = comb_m3_p2 / total_mag;     % [-3/+2]
    curve_data(7,k) = comb_m3_p3 / total_mag;     % [-3/+3]

end

% -------------------------------------------------------------------------
% 5. Plotting – match the structure of Peters' figure
% -------------------------------------------------------------------------

figure('Color','w','Position',[150 150 850 480]);
hold on; grid on; box on;

% Empirical split points where branches appear to merge/split in the figure
split1 = 0.40;   % where [-1] and [0] become [-1/+0]
split2 = 1.05;   % where [-2/+2] becomes important

m1 = mu_vals < split1;
m2 = (mu_vals >= split1) & (mu_vals < split2);
m3 = mu_vals >= split2;

% Segment 1: low mu – separate [-1] and [0]
plot(mu_vals(m1), curve_data(1,m1), 'k', 'LineWidth', 1.6);  % [-1]
plot(mu_vals(m1), curve_data(2,m1), 'k', 'LineWidth', 1.6);  % [+0]

% Segment 2: mid mu – combined [-1/+0]
plot(mu_vals(m2), curve_data(3,m2), 'k', 'LineWidth', 1.6);  % [-1/+0]

% Segment 3: higher mu – primarily [-2/+2] "fan" plus others
plot(mu_vals(m3), curve_data(5,m3), 'k', 'LineWidth', 1.6);  % [-2/+2]

% Higher-order combinations that appear later in Peters' figure
plot(mu_vals(mu_vals > split1), curve_data(4, mu_vals > split1), 'k', 'LineWidth', 1.4); % [-2/+1]
plot(mu_vals(mu_vals > 1.2),   curve_data(6, mu_vals > 1.2),   'k', 'LineWidth', 1.4); % [-3/+2]
plot(mu_vals(mu_vals > 1.8),   curve_data(7, mu_vals > 1.8),   'k', 'LineWidth', 1.4); % [-3/+3]

% Axes and labels
axis([0 mu_end 0 0.5]);
set(gca, 'FontSize', 12, 'LineWidth', 1.2, 'TickLabelInterpreter', 'latex');
xlabel('$\mu$', 'Interpreter','latex','FontSize',16);
ylabel('Modal Participation', 'Interpreter','latex','FontSize',14);
title('Harmonic Participation, Rotor-Blade Flapping ($p=1.0, \gamma=12$)', ...
      'Interpreter','latex','FontSize',14);

% Manual annotations (adjust positions to match your exact scan)
text(0.18, 0.42, '[-1]',      'Interpreter','latex','FontSize',12);
text(0.55, 0.28, '[+0]',      'Interpreter','latex','FontSize',12);
text(1.20, 0.35, '[-1/+0]',   'Interpreter','latex','FontSize',12);
text(1.10, 0.20, '[-2/+1]',   'Interpreter','latex','FontSize',12);
text(2.40, 0.16, '[-2/+2]',   'Interpreter','latex','FontSize',12);
text(2.40, 0.10, '[-3/+2]',   'Interpreter','latex','FontSize',12);
text(2.40, 0.06, '[-3/+3]',   'Interpreter','latex','FontSize',12);

hold off;
