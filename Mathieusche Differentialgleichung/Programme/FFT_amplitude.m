% Visualization of FFT Amplitudes for w = 0.3 (Perfect Visual & Scaling Corrections)
clear; clc; close all;

% --- Parameters ---
w = 0.3;            
Omega = 1.0;        
T = 2*pi / Omega;   
eps_vals = [0, 0.15, 0.5, 3.25, 3.8];
N_FFT = 4096;       

% --- Setup for Figure Saving ---
fDir = 'figureFolder'; % Folder for figures
if ~isdir(fDir) %#ok<*ISDIR>
    mkdir(fDir)
end
fDirPeters = fullfile(fDir,'FFT_Amplitude_SVG'); % Subfolder specific to Peters' plots
if ~isdir(fDirPeters)
    mkdir(fDirPeters)
end

% Setup identical color mapping matrix from the harmonic participation code (m_range = -3:3)
m_range = -3:3;
colors = lines(length(m_range)); 

figure('Color', 'w', 'Units', 'pixels', 'Position', [100, 100, 900, 950]);

for i = 1:length(eps_vals)
    epsilon = eps_vals(i);
    
    % --- 1. Solve for the State-Transition Matrix Phi(t) with strict tolerances ---
    options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
    D_func = @(t) [0, 1; -(w^2 + epsilon*sin(Omega*t)), 0];
    sol = ode45(@(t, x) reshape(D_func(t) * reshape(x, 2, 2), 4, 1), [0, T], reshape(eye(2), 4, 1), options);
    Phi_T = reshape(sol.y(:, end), 2, 2);
    
    % --- 2. Eigenanalysis ---
    [V, Lambda_mat] = eig(Phi_T);
    Lambda = diag(Lambda_mat);
    eta = log(Lambda) / T;
    [~, mode_idx] = max(imag(eta));
    eta_mode = eta(mode_idx);
    v_mode = V(:, mode_idx);
    
    % --- 3. Periodic Part A(t) and FFT ---
    t_fft = linspace(0, T, N_FFT + 1);
    t_fft(end) = []; 
    Phi_t_sampled = deval(sol, t_fft);
    
    Q_t_disp = zeros(1, N_FFT);
    for j = 1:N_FFT
        Phi_t = reshape(Phi_t_sampled(:, j), 2, 2);
        A_t = Phi_t * v_mode * exp(-eta_mode * t_fft(j));
        Q_t_disp(j) = A_t(1); 
    end
    
    % --- 4. FFT Execution & SCALING ---
    C_raw = fftshift(fft(Q_t_disp) / N_FFT);
    mag_C_raw = abs(C_raw);
    
    % Normalize coefficients by total energy so peaks never exceed 1.0 (Fix for epsilon = 3.25)
    total_energy = sum(mag_C_raw);
    if total_energy > 1e-12
        mag_C = mag_C_raw / total_energy; 
    else
        mag_C = mag_C_raw;
    end
    
    m_axis = (-N_FFT/2 : N_FFT/2 - 1);
    
    % --- 5. Plotting ---
    subplot(length(eps_vals), 1, i);
    hold on;
    
    % Main line is kept clean, uniform neutral gray
    plot(m_axis, mag_C, 'LineWidth', 1.5, 'Color', [0.7 0.7 0.7]);
    
    % Find the location of valid signal peaks
    [pks, locs] = findpeaks(mag_C, m_axis, 'MinPeakHeight', 0.015);
    
    % Only the peak markers are isolated and color-matched to plot A
    for k = 1:length(pks)
        m_val = round(locs(k)); % Round peak index location to its discrete integer m
        
        % Map the integer value back to the color matrix row index
        color_idx = find(m_range == m_val);
        if ~isempty(color_idx)
            peak_color = colors(color_idx, :); % Assign exact color matching participation curve
        else
            peak_color = [0.3 0.3 0.3]; % Dark gray fallback for peaks outside the tracking range
        end
        
        % Plot the isolated peak marker using its specific harmonic color tracking
        plot(locs(k), pks(k), 'o', 'MarkerFaceColor', peak_color, 'MarkerEdgeColor', 'k', 'MarkerSize', 6.5, 'LineWidth', 1);
        
        % Text labels for strong peaks
        if pks(k) > 0.04
            text(locs(k), pks(k) + 0.08, sprintf('$m=%d$', m_val), ...
                'HorizontalAlignment', 'center', 'FontSize', 9, 'Interpreter', 'latex');
        end
    end
    
    % Dynamic title text with parameter updates
    title(sprintf('$\\epsilon = %g, \\quad \\omega = %g$', epsilon, w), 'Interpreter', 'latex', 'FontSize', 12);
    
    ylabel('$|\phi_m|$', 'Interpreter', 'latex', 'FontSize', 11);
    grid on;
    xlim([-5, 5]);
    ylim([0, 1.1]); % Scaled limit down to 1.1 max boundary
    
    % Keep crucial tick evaluations (0, -1) cleanly visible layout-wise
    set(gca, 'XTick', -5:1:5, 'TickLabelInterpreter', 'latex');
    if i < length(eps_vals)
        set(gca, 'XTickLabel', []); % Clear intermediate clutter labels while keeping structural tick vectors
    end
end

xlabel('Harmonic Index ($m$)', 'Interpreter', 'latex', 'FontSize', 12);
sgtitle(sprintf('Harmonic Participation Spectrum for $\\omega = %g$', w), 'Interpreter', 'latex', 'FontSize', 14);

% Filename generation for saving the plot
pngname = 'FFT_Amplitude';
svgfile = fullfile(fDirPeters,[pngname,'.svg']);
print(svgfile, '-dsvg');