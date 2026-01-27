clear; clc; close all;

%% Global Parameters
Omega = 1;           
T     = 2*pi/Omega;  
w_values = 0.3; 
N_FFT  = 4096;
m_range = -5:5;      
x0 = eye(2);         
D = 0; % Damping (Set to 0 as per your previous prompt)

% ODE Matrix: x'' + 2Dx' + (eps + eps*cos(Om*t))x = 0
% Note: Using epsilon for both static and time-varying parts as per your snippet
ode_mat = @(t, nu0) [0, 1; -(nu0 + nu0*cos(Omega*t)), -2*D];

for w = w_values
    eps_end = 5; 
    eps_vals = linspace(0, eps_end, 500); 
    
    % Storage
    participation_data = zeros(length(eps_vals), length(m_range));
    composite_freq = zeros(length(eps_vals), 1); 
    growth_rate = zeros(length(eps_vals), 1); % Real part of exponent
    
    % Initial Guess: Start with the natural frequency
    % For the Mathieu equation, initial imag part is sqrt(static_term)
    eta_prev = 1i * sqrt(eps_vals(1)); 

    for k = 1:length(eps_vals)
      nu0 = eps_vals(k);
        
        % 1. Solve for Monodromy Matrix
        sol = ode45(@(t, x) reshape(ode_mat(t, nu0)*reshape(x,2,2),4,1), [0, T], reshape(x0,4,1));
        Phi_T = reshape(deval(sol, T), 2, 2);
        
        % 2. Extract Floquet Exponents
        [V, L_mat] = eig(Phi_T);
        eta_vals = log(diag(L_mat)) / T;
        
        % Branch Tracking: Choose the eigenvalue closest to the previous step
        [~, idx] = min(abs(eta_vals - eta_prev)); 
        eta_mode = eta_vals(idx);
        v_mode   = V(:, idx);
        eta_prev = eta_mode; 
        
        growth_rate(k) = real(eta_mode);
        
        % 3. Extract Periodic Part P(t) for Harmonic Analysis
        t_vec = linspace(0, T, N_FFT+1); t_vec(end) = [];
        Phi_t_all = deval(sol, t_vec);
        Q_t = zeros(N_FFT, 1);
        for j = 1:N_FFT
            Phi_curr = reshape(Phi_t_all(:,j), 2, 2);
            % P(t) = Phi(t) * v * exp(-eta*t)
            A_t = Phi_curr * v_mode * exp(-eta_mode * t_vec(j));
            Q_t(j) = A_t(1); 
        end
        
        C = fftshift(fft(Q_t)/N_FFT);
        freqs_fft_m = (-(N_FFT/2) : (N_FFT/2-1)); 
        
        % 4. Weighted Frequency Calculation
        % Normalize imag part to [-0.5, 0.5] range relative to Omega
        basis_imag = imag(eta_mode) / Omega;
        m_mags = zeros(1, length(m_range));
        branch_freqs = zeros(1, length(m_range));
        
        for i_m = 1:length(m_range)
            m = m_range(i_m);
            [~, f_idx] = min(abs(freqs_fft_m - m));
            m_mags(i_m) = abs(C(f_idx));
            % The physical frequency of the m-th harmonic is |m*Omega + imag(eta)|
            branch_freqs(i_m) = abs(m + basis_imag);
        end
        
        % Calculate participation weights
        if sum(m_mags) > 1e-12
            weights = m_mags / sum(m_mags);
            participation_data(k, :) = weights;
            % The composite frequency is the weighted average of harmonic frequencies
            composite_freq(k) = sum(branch_freqs .* weights);
        end
    end
    
    %% --- Refined Plotting ---
    figure('Name', 'Floquet Analysis Results', 'Color', 'w', 'Position', [100, 100, 900, 700]);
    
    % Subplot 1: Composite Frequency & Growth Rate
    subplot(3,1,1);
    yyaxis left
    plot(eps_vals, composite_freq, 'k-', 'LineWidth', 2);
    ylabel('Freq \omega_{comp}/\Omega');
    hold on;
    yyaxis right
    plot(eps_vals, growth_rate, 'r--', 'LineWidth', 1.5);
    ylabel('Growth Rate Re(\eta)');
    title(['Frequency and Stability (w = ', num2str(w), ')']);
    grid on;

    % Subplot 2: Participation Weights (Smoother display)
    subplot(3,1,2);
    area(eps_vals, participation_data, 'LineStyle', 'none');
    title('Harmonic Participation (Area Map)');
    ylabel('Weight');
    legend(cellstr(num2str(m_range', 'm=%d')), 'Location', 'eastoutside');
    grid on;

    % Subplot 3: Participation (Lines)
    subplot(3,1,3); 
    plot(eps_vals, participation_data, 'LineWidth', 1);
    title('Harmonic Participation (Lines)');
    xlabel('\nu'); ylabel('Weight');
    grid on;
end