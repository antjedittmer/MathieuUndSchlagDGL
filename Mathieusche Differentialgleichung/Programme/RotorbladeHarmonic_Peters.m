% -------------------------------------------------------------------------
% Exact Harmonic Participation: Rotor Blade Flapping
% Replicating Figure 3 from Peters (JAHS 2011)
% -------------------------------------------------------------------------
clear; clc; close all;

% --- Simulation Parameters ---
p = 1.0;            % Natural flap frequency
gamma = 12.0;       % Lock number
N_mu = 300;         % Resolution for smooth bends
mu_end = 3.0;
mu_vals = linspace(0.001, mu_end, N_mu);

% Storage for curves
% 1: [-1], 2: [0], 3: [-1/+0], 4: [-2/+1], 5: [-2/+2], 6: [-3/+2], 7: [-3/+3]
curve_data = zeros(7, N_mu);

clear; clc; close all;

% --- Parameters from Peters (2011) ---
p = 1.0;            % Natural flap frequency
gamma = 12.0;       % Lock number
N_mu = 200;
mu_vals = linspace(0, 2.5, N_mu); 

% Storage for the specific branches shown in Figure 3
results = struct('m1',[], 'p0',[], 'comb10',[], 'comb21',[], 'comb22',[]);

for k = 1:N_mu
    mu = mu_vals(k);
    % Exact Flap Equation Matrix from Peters (Eq 2)
    % beta'' + C(t)beta' + K(t)beta = 0
    C_t = @(t) (gamma/8) * (1 + (4/3)*mu*sin(t));
    K_t = @(t) p^2 + (gamma/8) * ((4/3)*mu*cos(t) + mu^2*sin(2*t));
    
    A_func = @(t) [0, 1; -K_t(t), -C_t(t)];

    % 1. Compute Monodromy Matrix (Transition matrix over 2*pi)
    [~, X] = ode45(@(t,x) reshape(A_func(t)*reshape(x,2,2),4,1), [0 2*pi], reshape(eye(2),4,1));
    Phi_T = reshape(X(end,:), 2, 2);
    
    % 2. Eigen-analysis (Choose the flapping mode)
    [V, L] = eig(Phi_T);
    [~, idx] = max(abs(diag(L))); 
    v_mode = V(:, idx);
    lambda = L(idx, idx);
    eta = log(lambda)/(2*pi); % Floquet exponent
    
    % 3. Extract periodic part Q(t) over a fine grid
    n_fft = 512;
    t_span = linspace(0, 2*pi, n_fft+1); t_span(end) = [];
    [~, X_path] = ode45(@(t,x) reshape(A_func(t)*reshape(x,2,2),4,1), t_span, reshape(eye(2),4,1));
    
    Q_t = zeros(n_fft, 1);
    for j = 1:n_fft
        Phi_t = reshape(X_path(j,:), 2, 2);
        % Floquet Theorem: x(t) = Q(t) * exp(eta*t)
        % Therefore: Q(t) = exp(-eta*t) * Phi(t) * v_mode
        temp_q = exp(-eta*t_span(j)) * (Phi_t * v_mode);
        Q_t(j) = temp_q(1);
    end
    
    % 4. FFT and Harmonic Extraction
    C = fftshift(fft(Q_t)/n_fft);
    mag = abs(C);
    mid = (n_fft/2) + 1; % Index of 0-th harmonic
    
    % Utility: Extract magnitude of harmonic 'm'
    h = @(m) mag(mid + m);
    total_energy = sum(mag); % Used for normalization
    
    % 5. Map the curves to the labels in the paper
    curve_data(1,k) = h(-1) / total_energy;         % [-1]
    curve_data(2,k) = h(0) / total_energy;          % [+0]
    curve_data(3,k) = (h(-1) + h(0)) / total_energy;% [-1/+0] (The combined bend)
    curve_data(4,k) = (h(-2) + h(1)) / total_energy;% [-2/+1]
    curve_data(5,k) = (h(-2) + h(2)) / total_energy;% [-2/+2]
    curve_data(6,k) = (h(-3) + h(2)) / total_energy;% [-3/+2]
    curve_data(7,k) = (h(-3) + h(3)) / total_energy;% [-3/+3]
end

% ================== PLOTTING ==================
figure('Color','w','Position',[150 150 850 480]);
hold on; grid on; box on;

% Transition Thresholds (Where the branches merge in the paper)
split1 = 0.40;  % Merging of -1 and 0
split2 = 1.05;  % Transition to -2/+2

% Plot segment 1: Low Advance Ratio (Separate -1 and 0)
m1 = mu_vals < split1;
plot(mu_vals(m1), curve_data(1,m1), 'k', 'LineWidth', 1.6);
plot(mu_vals(m1), curve_data(2,m1), 'k', 'LineWidth', 1.6);

% Plot segment 2: Mid Advance Ratio (Combined -1/+0)
m2 = (mu_vals >= split1) & (mu_vals < split2);
plot(mu_vals(m2), curve_data(3,m2), 'k', 'LineWidth', 1.6);

% Plot segment 3: High Advance Ratio (The "fanning out" of higher harmonics)
m3 = mu_vals >= split2;
plot(mu_vals(m3), curve_data(5,m3), 'k', 'LineWidth', 1.6);

% Higher order harmonics that appear early
plot(mu_vals(mu_vals > 0.4), curve_data(4, mu_vals > 0.4), 'k', 'LineWidth', 1.6);
plot(mu_vals(mu_vals > 1.2), curve_data(6, mu_vals > 1.2), 'k', 'LineWidth', 1.6);
plot(mu_vals(mu_vals > 1.8), curve_data(7, mu_vals > 1.8), 'k', 'LineWidth', 1.6);

% --- Styling and Labels ---
axis([0 3 0 0.5]);
set(gca, 'FontSize', 12, 'LineWidth', 1.2, 'TickLabelInterpreter', 'latex');
xlabel('$\mu$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('Modal Participation', 'Interpreter', 'latex', 'FontSize', 14);
title('Harmonic Participation (Peters 2011 Replication)', 'FontSize', 14);

% Manual Annotations
text(0.2, 0.42, '[-1]', 'FontSize', 12);
text(0.6, 0.28, '[+0]', 'FontSize', 12);
text(1.3, 0.35, '[-1/+0]', 'FontSize', 12);
text(1.2, 0.20, '[-2/+1]', 'FontSize', 12);
text(2.5, 0.15, '[-2/+2]', 'FontSize', 12);

    % Corrected periodic coefficients for Rotor Flapping
    C_t = @(t) (gamma/8) * (1 + (4/3)*mu*sin(t));
    K_t = @(t) p^2 + (gamma/8) * ((4/3)*mu*cos(t) + (mu^2)*sin(2*t)); 
    
    % State-space: x' = [0 1; -K -C]x
    A_func = @(t) [0, 1; -K_t(t), -C_t(t)];
    
    % 1. Compute Monodromy Matrix
    Phi_T = eye(2);
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
    
    % 4. Harmonic Analysis
    c = fftshift(fft(q_vals))/n_pts;
    freqs = (-n_pts/2 : n_pts/2 - 1);
    mags = abs(c);
    
    % Helper to get magnitude of harmonic 'n'
    get_mag = @(n) mags(freqs == n);
    
    % 5. Store data for plotting
    % Peters plots relative participation (normalized by sum of magnitudes)
    total = sum(mags);
    results.m1(k) = get_mag(-1)/total;
    results.p0(k) = get_mag(0)/total;
    results.comb10(k) = (get_mag(-1) + get_mag(0))/total;
    results.comb21(k) = (get_mag(-2) + get_mag(1))/total;
    results.comb22(k) = (get_mag(-2) + get_mag(2))/total;

% --- Plotting to match JAHS 2011 Figure 3 ---
figure('Color','w'); hold on; grid on;
% In the paper, curves are often solid. We use thresholds to show the 'splits'.
% The [-1] and [0] branches start separate then merge.
plot(mu_vals, results.m1, 'k', 'LineWidth', 1.5);
plot(mu_vals, results.p0, 'k', 'LineWidth', 1.5);
plot(mu_vals, results.comb21, 'k--', 'LineWidth', 1.2); % Higher order sidebands

% Styling
xlabel('Advance Ratio \mu');
ylabel('Harmonic Participation \phi_n');
title('Replication of Peters (2011) Fig. 10');
axis([0 2.5 0 0.6]);

