clear; clc; close all;

%% === SETUP ===
dDirA = fullfile('dataFolder', 'dataFolder_Arnold_Classic_Symmetric_test');

%% === 1. LOAD ARNOLD  ===
matA = 'STRUTTscheKarte_D1dot5e-01_SW1dot0e-01_unt0.mat';
S_A = load(fullfile(dDirA, matA));
nu_A = S_A.CharEx(:,1);
ImA1 = S_A.CharEx(:,7);  % Trusted Im(s_R1)
ImA2 = S_A.CharEx(:,8);  % Trusted Im(s_R2)
fprintf('Arnold: %d points ✓\n', length(nu_A));

%% Global Parameters
Omega = 1;           
T     = 2*pi/Omega;  
nu_vals = linspace(0.01, 9, 500); 
m_range = -5:5;                   
N_FFT  = 2048;
x0 = eye(2);         
D = 0;                            

% Pre-allocate storage
growth_rate        = zeros(length(nu_vals), 1); 
imag_exponent      = zeros(length(nu_vals), 1); 
participation_data = zeros(length(nu_vals), length(m_range));
branch_freqs_all   = zeros(length(nu_vals), length(m_range));
composite_freq     = zeros(length(nu_vals), 1);

% Initial Guess for Tracking
eta_prev = 1i * sqrt(nu_vals(1)); 

%% Computation Loop
for k = 1:length(nu_vals)
    nu = nu_vals(k);
    
    % 1. Solve for Monodromy Matrix
    ode_mat = @(t, x) [0, 1; -(nu + nu*cos(Omega*t)), -2*D] * reshape(x, 2, 2);
    [~, sol_raw] = ode45(@(t, x) reshape(ode_mat(t, x), 4, 1), [0, T], reshape(x0, 4, 1));
    Phi_T = reshape(sol_raw(end, :), 2, 2);
    
    % 2. Extract Floquet Exponents
    [V, L_mat] = eig(Phi_T);
    eta_vals = log(diag(L_mat)) / T;
    
    [~, idx] = min(abs(eta_vals - eta_prev)); 
    eta_mode = eta_vals(idx);
    v_mode   = V(:, idx);
    eta_prev = eta_mode; 
    
    growth_rate(k)   = real(eta_mode);
    imag_exponent(k) = imag(eta_mode);
    
    % 3. Extract Periodic Part P(t) via FFT
    t_fft = linspace(0, T, N_FFT+1); t_fft(end) = [];
    sol_obj = ode45(@(t, x) reshape(ode_mat(t, x), 4, 1), [0, T], reshape(x0, 4, 1));
    Phi_t_steps = deval(sol_obj, t_fft);
    
    Q_t = zeros(N_FFT, 1);
    for j = 1:N_FFT
        Phi_curr = reshape(Phi_t_steps(:,j), 2, 2);
        P_t = Phi_curr * v_mode * exp(-eta_mode * t_fft(j));
        Q_t(j) = P_t(1); 
    end
    
    C = fftshift(fft(Q_t)/N_FFT);
    freqs_fft_m = (-(N_FFT/2) : (N_FFT/2-1)); 
    
    % 4. Frequency and Participation
    for i_m = 1:length(m_range)
        m = m_range(i_m);
        [~, f_idx] = min(abs(freqs_fft_m - m));
        participation_data(k, i_m) = abs(C(f_idx));
        branch_freqs_all(k, i_m) = abs(m*Omega + imag(eta_mode));
    end
    
    % Normalize Weights and calculate composite freq
    weights = participation_data(k, :);
    total_mag = sum(weights);
    if total_mag > 1e-12
        weights = weights / total_mag;
        participation_data(k, :) = weights;
        composite_freq(k) = sum(branch_freqs_all(k, :) .* weights);
    end
end

%% --- Refined Plotting ---
figure('Name', 'Floquet Analysis (nu 0-9)', 'Color', 'w', 'Position', [100, 50, 900, 950]);

% Subplot 1: Composite Frequency & Growth Rate
subplot(4,1,1);
yyaxis left
plot(nu_vals, composite_freq, 'k-', 'LineWidth', 2);
ylabel('Freq \omega_{comp}/\Omega');
hold on;
plot(nu_A,ImA2)

yyaxis right
plot(nu_vals, growth_rate, 'r--', 'LineWidth', 1.5);
ylabel('Growth Rate Re(\eta)');
title('Frequency and Stability');
grid on;

% Plot 2: Imaginary Part
subplot(4,1,2);
ImA2Interp = interp1(nu_A,ImA2,nu_vals);
plot(nu_vals, composite_freq -ImA2Interp', 'Color', [0, 0.5, 0], 'LineWidth', 1.5);
ylabel('Im(\eta)');
title('Subplot 2: Imaginary Part of Characteristic Exponent');
grid on;

% Plot 3: Frequency Branches
subplot(4,1,3);
plot(nu_vals, branch_freqs_all, 'LineWidth', 1);
ylabel('Branch Freqs');
title('Subplot 3: Individual Frequency Branches');
grid on;

% Plot 4: Harmonic Participation
subplot(4,1,4);
set(gca, 'ColorOrder', jet(length(m_range))); 
plot(nu_vals, participation_data, 'LineWidth', 1.2);
xlabel('\nu (Addition Factor)');
ylabel('Weight');
title('Subplot 4: Harmonic Participation');
grid on;