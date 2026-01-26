% Refactored Script: Floquet Analysis of Mathieu's Equation
% Including Weighted Composite Frequency Calculation
clear; clc; close all;

%% ------------------------------------------------------------------------
% Global Parameters
%% ------------------------------------------------------------------------
Omega = 1;           
T     = 2*pi/Omega;  
w_values = 0.3 %  [0.3, 0.5, 0.7];
N_FFT  = 4096;
m_range = -5:5;      
 
x0 = eye(2);         
D = 0
ode_mat = @(t, epsilon, w_sq) [0, 1; -(w_sq + epsilon*sin(Omega*t)), -2*D];

%% ========================================================================
% Loop over w-values
%% ========================================================================
for w = w_values
    w_sq = w^2;
    
    % Setup Epsilon Range
    eps_end = 3.5; if abs(w - 0.3) < 1e-3, eps_end = 5; end
    eps_vals = linspace(0, eps_end, 1000);
    
    % Pre-allocate storage
    freq_results = struct();
    for m = m_range
        f_field = matlab.lang.makeValidName(['m_' num2str(m)]);
        freq_results.(f_field) = zeros(length(eps_vals), 2);
    end
    participation_data = zeros(length(eps_vals), length(m_range));
    composite_freq = zeros(length(eps_vals), 1); % Store weighted sum
    
    % Mode tracking initialization
    eta_prev = -w * Omega;
    fprintf('Processing w = %.1f...\n', w);
    
    for k = 1:length(eps_vals)
        epsilon = eps_vals(k);
        
        % 1. Solve ODE once
        sol = ode45(@(t, x) reshape(ode_mat(t, epsilon, w_sq)*reshape(x,2,2),4,1), [0, T], reshape(x0,4,1));
        Phi_T = reshape(deval(sol, T), 2, 2);
        
        % 2. Floquet Exponents
        [V, L_mat] = eig(Phi_T);
        eta_vals = log(diag(L_mat)) / T;
        [~, idx] = min(abs(eta_vals - eta_prev)); 
        eta_mode = eta_vals(idx);
        v_mode   = V(:, idx);
        eta_prev = eta_mode; 
        
        % 3. Harmonic Content (FFT)
        t_vec = linspace(0, T, N_FFT+1); t_vec(end) = [];
        Phi_t_all = deval(sol, t_vec);
        Q_t = zeros(N_FFT, 1);
        for j = 1:N_FFT
            Phi_curr = reshape(Phi_t_all(:,j), 2, 2);
            A_t = Phi_curr * v_mode * exp(-eta_mode * t_vec(j));
            Q_t(j) = A_t(1); 
        end
        
        C = fftshift(fft(Q_t)/N_FFT);
        freqs_fft = (-(N_FFT/2) : (N_FFT/2-1)) / T;
        
        % 4. Calculate Individual Branch Frequencies and Modal Participation
        basis_imag = imag(eta_mode) / Omega;
        current_step_branch_freqs = zeros(1, length(m_range));
        m_mags = zeros(1, length(m_range));
        i_m = 0;
        for m = m_range
            i_m = i_m + 1;
            f_field = matlab.lang.makeValidName(['m_' num2str(m)]);
            branch_f = (m == 0) * basis_imag + (m ~= 0) * (abs(m)*Omega + sign(m)*basis_imag);
            freq_results.(f_field)(k, :) = [epsilon, branch_f];
            [~, f_idx] = min(abs(freqs_fft - m/T));
            m_mags(i_m) = abs(C(f_idx));

            % Identify the specific frequency for this harmonic branch at current epsilon
            current_step_branch_freqs(i_m) = (m == 0) * basis_imag + ...
                (m ~= 0) * (abs(m)*Omega + sign(m)*basis_imag);
        end
        
        if sum(m_mags) > 1e-12
            weights = m_mags / sum(m_mags);
            participation_data(k, :) = weights;
            
            % THE MULTIPLICATION AND ADDITION STEP:
            % Composite Frequency = Sum( Branch_Frequency_i * Participation_Weight_i )
            composite_freq(k) = sum(current_step_branch_freqs .* weights);
        end
    end
    
    %% --- Plotting Results ---
    
    % 1. Composite Frequency Plot (The new requested calculation)
    figure('Name', sprintf('Composite Freq w=%.1f', w), 'Color', 'w');
    plot(eps_vals, composite_freq, 'k-', 'LineWidth', 2);
    title(['Weighted Composite Frequency vs \epsilon (w = ', num2str(w), ')']);
    xlabel('\epsilon'); ylabel('Weighted \omega / \Omega'); grid on;

    % 2. Participation Plot
    figure('Name', sprintf('Participation w=%.1f', w), 'Color', 'w');
    plot(eps_vals, participation_data, 'LineWidth', 1.5);
    title(['Harmonic Participation for w = ', num2str(w)]);
    xlabel('\epsilon'); ylabel('Participation'); grid on;
    legend(cellstr(num2str(m_range', 'm=%d')), 'Location', 'bestoutside');

    % 3. Frequency Branches Plot
    figure('Name', sprintf('Freq Branches w=%.1f', w), 'Color', 'w'); hold on;
    f_names = fieldnames(freq_results);
    for i = 1:length(f_names)
        data = freq_results.(f_names{i});
        plot(data(:,1), data(:,2), '.', 'MarkerSize', 4);
    end
    title(['Individual Frequency Branches for w = ', num2str(w)]);
    xlabel('\epsilon'); ylabel('\omega / \Omega'); grid on;
end