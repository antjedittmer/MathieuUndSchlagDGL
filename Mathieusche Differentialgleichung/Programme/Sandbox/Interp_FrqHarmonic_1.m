
% Floquet Analysis of Mathieu's Equation
% Consolidation of Frequency, Participation, and Waveform analysis.

clear; clc; close all;

%% ------------------------------------------------------------------------
% Global Parameters
%% ------------------------------------------------------------------------
Omega = 1;           
T     = 2*pi/Omega;  
w_values = [0.3, 0.5, 0.7];
N_FFT  = 4096;
m_range = -4:4;      % Combined range for frequency and participation
m_range_harm = -3:3; % Subset for modal participation plots
x0 = eye(2);         % Initial condition

% Single ODE definition: x'' + (w^2 + epsilon*sin(Omega t)) x = 0
ode_mat = @(t, epsilon, w_sq) [0, 1; -(w_sq + epsilon*sin(Omega*t)), 0];

%% ========================================================================
% Loop over w-values
%% ========================================================================
for w = w_values
    w_sq = w^2;
    
    % Setup Epsilon Range
    eps_end = 3.5; if abs(w - 0.3) < 1e-3, eps_end = 5; end
    eps_vals = linspace(0, eps_end, 400);
    
    % Pre-allocate storage
    freq_results = struct();
    for m = m_range
        f_field = matlab.lang.makeValidName(['m_' num2str(m)]);
        freq_results.(f_field) = zeros(length(eps_vals), 2);
    end
    participation_data = zeros(length(eps_vals), length(m_range_harm));
    
    % Mode tracking initialization
    eta_prev = -w * Omega;

    fprintf('Processing w = %.1f...\n', w);

    for k = 1:length(eps_vals)
        epsilon = eps_vals(k);

        % 1. Solve ODE once for this epsilon
        sol = ode45(@(t, x) reshape(ode_mat(t, epsilon, w_sq)*reshape(x,2,2),4,1), [0, T], reshape(x0,4,1));
        Phi_T = reshape(deval(sol, T), 2, 2);

        % 2. Floquet Exponents and Mode Tracking
        [V, L_mat] = eig(Phi_T);
        eta_vals = log(diag(L_mat)) / T;
        [~, idx] = min(abs(eta_vals - eta_prev)); % Track consistent branch
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
            Q_t(j) = A_t(1); % Displacement component
        end
        
        C = fftshift(fft(Q_t)/N_FFT);
        freqs = (-(N_FFT/2) : (N_FFT/2-1)) / T;

        % 4. Store Branch Frequencies
        basis_imag = imag(eta_mode) / Omega;
        for m = m_range
            f_field = matlab.lang.makeValidName(['m_' num2str(m)]);
            branch_f = (m == 0) * basis_imag + (m ~= 0) * (abs(m)*Omega + sign(m)*basis_imag);
            freq_results.(f_field)(k, :) = [epsilon, branch_f];
        end

        % 5. Store Participation
        m_mags = zeros(1, length(m_range_harm));
        for i_m = 1:length(m_range_harm)
            [~, f_idx] = min(abs(freqs - m_range_harm(i_m)/T));
            m_mags(i_m) = abs(C(f_idx));
        end
        if sum(m_mags) > 1e-12
            participation_data(k, :) = m_mags / sum(m_mags);
        end
        
        % 6. Sample Waveform (at epsilon approx 3.0)
        if k == round(length(eps_vals)*0.75) 
            figure('Name', sprintf('Waveform w=%.1f', w));
            plot(t_vec, real(Q_t), 'LineWidth', 1.5);
            title(['Waveform Q(t) at $\epsilon \approx$ ', num2str(epsilon)]);
            grid on; xlabel('Time'); ylabel('Q(t)');
        end
    end

    %% --- Plotting Results ---
    % Frequency Plot
    figure('Name', sprintf('Freq w=%.1f', w)); hold on;
    f_names = fieldnames(freq_results);
    for i = 1:length(f_names)
        data = freq_results.(f_names{i});
        plot(data(:,1), data(:,2), '.', 'MarkerSize', 4);
    end
    title(['Frequency vs \epsilon for w = ', num2str(w)]);
    xlabel('\epsilon'); ylabel('\omega / \Omega'); grid on;

    % Participation Plot
    figure('Name', sprintf('Participation w=%.1f', w));
    plot(eps_vals, participation_data, 'LineWidth', 1.5);
    title(['Harmonic Participation for w = ', num2str(w)]);
    xlabel('\epsilon'); ylabel('Participation'); grid on;
    legend(cellstr(num2str(m_range_harm', 'm=%d')));
end

