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
end

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