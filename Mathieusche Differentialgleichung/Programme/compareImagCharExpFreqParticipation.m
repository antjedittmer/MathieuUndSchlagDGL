clear; clc; close all;

%% === SETUP ===
dDirA = fullfile('dataFolder', 'dataFolder_Arnold_Classic_Symmetric_test');

%% === 1. LOAD ARNOLD REFERENCE DATA ===
matA = 'STRUTTscheKarte_D1dot5e-01_SW1dot0e-01_unt0.mat';
S_A = load(fullfile(dDirA, matA));
nu_A = S_A.CharEx(:,1);
ImA2 = S_A.CharEx(:,8);  % Trusted Im(s_R2) - Arnold reference
fprintf('Arnold: %d points ✓\n', length(nu_A));

%% === FIGURE FOLDER SETUP ===
fDir = 'figureFolder';
if ~isdir(fDir)
    mkdir(fDir);
end
fDirPeters = fullfile(fDir, 'figureFolderPeters');
if ~isdir(fDirPeters)
    mkdir(fDirPeters);
end

%% === GLOBAL PARAMETERS ===
Omega = 1;           
T = 2*pi/Omega;  
nu_vals = linspace(0.01, 9, 500); 
m_range = -4:4;                   
N_FFT = 2048;
D = 0;                            
x0 = eye(2);

%% === PRE-ALLOCATION ===
growth_rate = zeros(length(nu_vals), 1); 
composite_freq = zeros(length(nu_vals), 1); 
participation_data = zeros(length(nu_vals), length(m_range));
branch_freqs_all = zeros(length(nu_vals), length(m_range));

%% === INITIAL GUESS FOR MODE TRACKING ===
eta_prev = 1i * sqrt(nu_vals(1)); 

%% === MAIN COMPUTATION LOOP ===
for k = 1:length(nu_vals)
    nu = nu_vals(k);
    
    % Solve for Monodromy Matrix
    ode_mat = @(t, x) [0, 1; -(nu + nu*cos(Omega*t)), -2*D] * reshape(x, 2, 2);
    [~, sol_raw] = ode45(@(t, x) reshape(ode_mat(t, x), 4, 1), [0, T], reshape(x0, 4, 1));
    Phi_T = reshape(sol_raw(end, :), 2, 2);
    
    % Extract Floquet Exponents (Peters method)
    [V, L_mat] = eig(Phi_T);
    eta_vals = log(diag(L_mat)) / T;
    [~, idx] = min(abs(eta_vals - eta_prev)); 
    eta_mode = eta_vals(idx);
    v_mode = V(:, idx);
    eta_prev = eta_mode; 
    
    growth_rate(k) = real(eta_mode);
    
    % Extract Periodic Part P(t) via FFT for participation
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
    
    % Frequency and Participation Calculation
    for i_m = 1:length(m_range)
        m = m_range(i_m);
        [~, f_idx] = min(abs(freqs_fft_m - m));
        participation_data(k, i_m) = abs(C(f_idx));
        branch_freqs_all(k, i_m) = abs(m*Omega + imag(eta_mode));
    end
    
    % Normalize weights and compute composite frequency
    weights = participation_data(k, :);
    total_mag = sum(weights);
    if total_mag > 1e-12
        weights = weights / total_mag;
        participation_data(k, :) = weights;
        composite_freq(k) = sum(branch_freqs_all(k, :) .* weights);
    end
end

%% === PLOTTING ===
nc = size(unique(lines,'rows'),1); % Numbers of different colors: nc = 7
lines0 = lines;
cl = [lines0(1:nc,:); 0*ones(1,3); 0.5*ones(1,3)]; % 7 line colors, black and grey

fig = figure('Name', 'Characteristic Exponent: Peters vs Arnold', 'Color', 'w'); 
pos0 = get(0,'defaultFigurePosition');
fig.Position = [pos0(1), pos0(2)- 0.3*pos0(4), pos0(3), 1.5*pos0(4)];
tiledlayout(4,1,'TileSpacing','tight');


% Subplot 1: Composite Frequency (Peters solid) vs Arnold (dashed) + Growth Rate
nexttile; %subplot(4,1,1);
yyaxis left
plot(nu_vals, composite_freq, '-', 'LineWidth', 2, 'Color', cl(1,:), 'DisplayName', 'Peters calculation');
hold on;
plot(nu_A, ImA2, '--', 'LineWidth', 2, 'Color', cl(1,:), 'DisplayName', 'Arnold calculation');
ylabel('\omega = Im(s_R)');
yyaxis right
plot(nu_vals, growth_rate, '-', 'LineWidth', 1.5, 'Color', cl(2,:), 'DisplayName', 'Real Part Re(s_R)');
ylabel('\sigma = Re(s_R)');
title('Mathieu ODE: x''''(t)+ 2Dx''(t) +(\nu_0^2 + \nu_c^2 cos(t))x(t) = 0, D = 0, \nu_0 = \nu_c','Real and Imaginary Part of Characteristic Exponent');
legend('Location','best');
grid on;

% Subplot 2: Difference (Peters - Arnold)
nexttile; %subplot(4,1,2);
ImA2Interp = interp1(nu_A, ImA2, nu_vals);
plot(nu_vals, composite_freq - ImA2Interp', 'Color', cl(1,:), 'LineWidth', 1.5, 'DisplayName', 'Peters - Arnold');
ylabel('Δω');
title('Difference Imaginary Part of Characteristic Exponent');
legend('Location','best');
grid on;

% Subplot 3: Frequency Branches
nexttile; %subplot(4,1,3);
for idx = 1: size(branch_freqs_all,2)
    plot(nu_vals, branch_freqs_all(:,idx), 'LineWidth', 1, 'color',cl(idx,:));
    hold on;
end
ylabel('Branch Freqs');
title('Individual Frequency Branches');
legend(arrayfun(@(m) sprintf('m=%d', m_range(m)), 1:length(m_range), 'UniformOutput', false), 'Location','northeastoutside');
grid on;

% Subplot 4: Harmonic Participation
nexttile; %subplot(4,1,4);
%set(gca, 'ColorOrder', jet(length(m_range))); 
for idx = 1: size(branch_freqs_all,2)
    plot(nu_vals, participation_data(:,idx), 'LineWidth', 1, 'color',cl(idx,:));
    hold on;
end
xlabel('Amplification factor \nu^2_c');
ylabel('Weight');
title('Harmonic Participation');
legend(arrayfun(@(m) sprintf('m=%d', m_range(m)), 1:length(m_range), 'UniformOutput', false), 'Location','northeastoutside');
grid on;

%% === SAVE FIGURE ===
pngname = sprintf('compare_char_exp_freq_participation_w1dot0.png');
pngfile = fullfile(fDirPeters, pngname);
saveas(fig, pngfile);
fprintf('Figure saved: %s\n', pngfile);
