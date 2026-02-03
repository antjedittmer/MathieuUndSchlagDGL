% MATLAB code to recreate Figure 21: Harmonic participation factors for rotor 
% blade flapping motion (p=1.0, γ=9.6) using Floquet theory approach from Peters (2009)
% Solves periodic blade flapping equation and computes harmonic participation factors

clear; clc; close all;

%% Parameters (from paper context)
p = 1.0;           % Flap frequency parameter (1/rev)
gamma = 9.6;       % Lock number
mu_range = 0:0.01:2.0;  % Advance ratio range (matches Fig 21 x-axis)
N_harm = 4;        % Number of harmonics to consider [+/-1,+/-2,+/-3,+/-4]
dt = 0.01;         % Time step for one rotor revolution (T=2*pi)
t = 0:dt:2*pi-dt;  % One period (azimuth ψ from 0 to 2π)
n_mu = length(mu_range);

%% Pre-allocate participation factors φ_n (Fig 21 y-axis data structure)
phi = zeros(n_mu, 2*N_harm+1);  % Columns: [+0], [+1],...,[+N], [-1],...,[-N]
labels = cell(1, size(phi,2));
labels{1} = '[+0]'; labels{end} = '[-4/+4]';  % Simplified labeling

%% Main loop: Compute Floquet solution for each advance ratio μ
for i = 1:n_mu
    mu = mu_range(i);
    
    % Periodic flap equation coefficients (simplified Peters formulation)
    % β'' + γ*(aerodynamic damping terms) + p^2 β = periodic forcing
    K_t = p^2 * ones(size(t));                    % Stiffness (flap freq^2)
    C_t = gamma * 0.5 * (1 + mu*cos(t));          % Damping (inflow variation)
    M_t = ones(size(t));                          % Mass = 1 (normalized)
    
    % Assemble state-space: [β; β'] → A(t)[x] = [x]'
    A = zeros(2,2,length(t));
    for k=1:length(t)
        A(:,:,k) = [-C_t(k)/M_t(k), -K_t(k)/M_t(k); 1, 0];
    end
    
    % Numerical transition matrix over one period using state transition
    Phi_T = eye(2);  % Floquet transition matrix Φ(T)
    n_steps = length(t);
    for k=1:n_steps-1
        Phi_step = expm(A(:,:,k)*dt);  % Local state transition
        Phi_T = Phi_step * Phi_T;
    end
    
    % Floquet exponents η = (1/T) log(eigen(Φ(T)))
    T = 2*pi;
    [V,D] = eig(Phi_T);
    eta = log(diag(D))/T;  % Complex exponents σ + iω
    
    % Select primary flap mode (closest to p=1.0, real part near 0)
    [~, idx] = min(abs(imag(eta) - p) + abs(real(eta)));
    eta_mode = eta(idx);
    
    % Periodic eigenvector A(t) via Eq(14): A(t) = Φ(t) * V(:,idx) * exp(-η t)
    A0 = V(:,idx);  % A(0)
    x_mode = zeros(2,length(t));  % Reconstruct full mode shape
    Phi_t = eye(2);
    for k=1:length(t)
        Phi_step = expm(A(:,:,k)*dt);
        Phi_t = Phi_step * Phi_t;
        x_mode(:,k) = Phi_t * A0 * exp(-eta_mode * t(k));
    end
    
    % Extract flapping β(t) = x_mode(1,:), normalize
    beta_t = x_mode(1,:);
    beta_t = beta_t / max(abs(beta_t));
    
    % Fourier analysis: Compute harmonic participation φ_n
    % Use FFT to get coefficients for n = -N_harm to +N_harm
    N_FFT = 2^nextpow2(length(beta_t));
    fft_beta = fft(beta_t, N_FFT);
    f = (0:N_FFT-1)/N_FFT * (1/dt);  % Frequency axis (rad/s, 1/rev normalized)
    
    % Participation factors |c_n|^2 / sum(|c_n|^2) for key harmonics
    harm_idx = [-N_harm:N_harm] + 1;  % FFT indices for harmonics
    c_n = fft_beta(harm_idx);
    phi_normsq = abs(c_n).^2;
    phi(i,1:length(harm_idx)) = phi_normsq / sum(phi_normsq);
end

%% Plot Figure 21 recreation
figure('Position',[100 100 1000 600]);
plot(mu_range, phi, 'LineWidth', 1.5);
xlabel('Advance Ratio \mu');
ylabel('Harmonic Participation Factor \phi_n');
title('Figure 21 Recreation: Rotor Blade Flapping (p=1.0, \gamma=9.6)');
legend(labels, 'Location','best');
grid on; xlim([0 2]); ylim([0 0.5]);
