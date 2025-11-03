% MATLAB Code to Recreate Mathieu's Equation Frequency Plot (Figure 2 Appearance)
% Separates branches according to the integer multiple 'm' used.
clear; clc; close all;

% --- Setup for Figure Saving ---
fDir = 'figureFolder'; % Folder for figures
if ~isdir(fDir) %#ok<ISDIR>
    mkdir(fDir)
end
fDirPeters = fullfile(fDir,'figureFolderPeters'); % Subfolder specific to Peters' plots
if ~isdir(fDirPeters) %#ok<ISDIR>
    mkdir(fDirPeters)
end

% Plotting style selection
K = 'ColoredLines';
% K = 'BlackLines';
useK = strcmp(K,'BlackLines');

% -----------------------------------------------------------------------
% --- Parameters and Initialization ---
% -----------------------------------------------------------------------
w_sq = 0.3^2;    % Natural frequency squared (w^2 in the ODE)
w = sqrt(w_sq);  % Natural frequency w = 0.7
Omega = 1;       % Fundamental frequency (Omega = 1 rad/s)
% Define base frequency omega0 based on Peters' convention (0 <= omega0 <= Omega/2).
% For w=0.7, the base frequency omega0 is 0.3.
omega0 = 0.3; 

% Filename generation for saving the plot
pngname = strrep(sprintf('PetersFrequency%s_%2.1f',K,w),'.','dot');
pngfile = fullfile(fDirPeters,[pngname,'.png']);

% Determine the range of epsilon (x-axis limit)
if abs(w - 0.3) < 0.001
    eps_end = 5; 
    eps_no = 1000; % High resolution for w=0.3
else
    eps_end = 3.5; % Axis limit for w=0.7 plot (as in the JAHS2011 paper)
    eps_no = 150; 
end

eps_vals = linspace(0, eps_end, eps_no);
m_range = (-4:4); % Integer multiple range for plotting branches (m*Omega)

% Structure to hold results, organized by branch 'm' (e.g., m_0, m_neg_1)
results_by_branch = struct();
for m = m_range
    % Create a valid field name (MATLAB field names cannot start with a sign)
    if m < 0
        field_name = ['m_neg_', num2str(abs(m))];
    else
        field_name = ['m_', num2str(m)];
    end
    results_by_branch.(field_name) = []; % Initialize with empty array
end

% -----------------------------------------------------------------------
% --- Floquet Exponent Calculation and Branch Separation ---
% -----------------------------------------------------------------------
T = 2*pi/Omega; % Period
x0 = eye(2); % Intial state matrix for Phi(0) = I

for k = 1:length(eps_vals)
    epsilon = eps_vals(k);
    
    % The state matrix D(t) for the ODE: ddot(x) + (w^2 + epsilon*sin(Omega*t)) * x = 0
    D_func = @(t) [0, 1; -(w_sq + epsilon*sin(t)), 0];
    
    % Solve for the Transition Matrix Phi(T)
    [~, Phi_t] = ode45(@(t, x) reshape(D_func(t) * reshape(x, 2, 2), 4, 1), [0, T], reshape(x0, 4, 1));
    Phi_T = reshape(Phi_t(end, :), 2, 2); % Monodromy Matrix at T
    
    % Calculate Floquet Exponents: eta = log(Lambda) / T
    Lambda = eig(Phi_T);
    eta = log(Lambda) / T;
    
    % Extract the Imaginary Part (normalized frequency, mu*Omega = Im(eta))
    % omega/Omega = Im(eta)/Omega (where Omega=1, so omega = Im(eta))
    normalized_omega = imag(eta) / Omega;

    % --- Separate and Store Branches (m) ---
    for r = 1:length(normalized_omega)
        freq_r = normalized_omega(r);
        
        % 1. Find the basis frequency (principal value in the range [-0.5, 0.5])
        % This step effectively extracts the fractional part of the frequency, 
        % which corresponds to the base frequency component (omega0/Omega).
        basis_freq = mod(freq_r + 0.5, 1) - 0.5;
        
        % 2. Map this base frequency to all possible branches 'm'
        % The true frequency can be shifted by any integer multiple 'm' of Omega.
        % Freq_m = basis_freq + m
        for m = m_range
            branch_freq = basis_freq + m;
            
            % Determine the valid field name
            if m < 0
                field_name = ['m_neg_', num2str(abs(m))];
            else
                field_name = ['m_', num2str(m)];
            end
            
            % Store [epsilon, frequency] for this specific branch
            results_by_branch.(field_name) = [results_by_branch.(field_name); epsilon, branch_freq];
        end
    end
end

% -----------------------------------------------------------------------
% --- Plotting ---
% -----------------------------------------------------------------------
figure;
hold on;
color_map = lines(length(m_range)); % Color map for distinct lines

% Title Update: Includes ODE and w for clarity
ode_str = '$\ddot{x} + (w^2 + \epsilon\sin(\Omega t)) x = 0$';
w_str = num2str(w, '%1.1f');
new_title = ['Frequency vs. $\epsilon$, ', ode_str, ', $w = ', w_str, '$'];
title(new_title, 'Interpreter', 'latex');

xlabel('$\epsilon$', 'FontSize', 14, 'Interpreter', 'latex');
% YLabel reflects normalized frequency
ylabel('Frequency ($\omega/\Omega$)', 'FontSize', 14, 'Interpreter', 'latex');

idx = 1;
for m = m_range
    % Determine field name
    if m < 0
        field_name = ['m_neg_', num2str(abs(m))];
    else
        field_name = ['m_', num2str(m)];
    end

    % --- Legend Calculation: Peters' Frequency Convention ---
    % m corresponds to the frequency shift from the base frequency omega0:
    % For m >= 0: omega/Omega approx (m + omega0/Omega)
    % For m < 0: omega/Omega approx (|m| - omega0/Omega)
    if m >= 0
        freq_normalized = omega0/Omega + m;
        freq_str = sprintf('$m=%+d \\rightarrow \\omega/\\Omega \\approx %.1f$', m, freq_normalized);
    elseif m < 0
        m_abs = abs(m);
        freq_normalized = m_abs - omega0/Omega;
        freq_str = sprintf('$m=%+d \\rightarrow \\omega/\\Omega \\approx %.1f$', m, freq_normalized);
    end

    % Get data for plotting
    data = results_by_branch.(field_name);
    
    % Plot using dots/scatters to form the continuous curves
    if useK == 1
        plot(data(:, 1), data(:, 2), '.', 'Color','k', 'MarkerSize', 8, 'DisplayName', freq_str);
    else
        plot(data(:, 1), data(:, 2), '.', 'Color', color_map(idx, :), 'MarkerSize', 8, 'DisplayName', freq_str);
    end
    idx = idx + 1;
end

% Set Axis limits
axis([0 eps_end 0 4.5]);
grid on;
set(gca, 'TickLabelInterpreter', 'latex');

% Add Legend if colors are used
if ~useK
    legend('Location', 'northeastoutside', 'Interpreter', 'latex');
end
hold off;

% -----------------------------------------------------------------------
% --- Branch Label Annotations ---
% -----------------------------------------------------------------------
% These labels mark the primary frequency branches at different regions
% of the stability chart.
hold on;
% Left-side labels (ε -> 0)
text(0.05, 0.05, '[+0]', 'FontSize', 10, 'Interpreter', 'latex'); % Near w/Omega = 0.7 - 0.5 = 0.2 (using basis freq)
text(0.05, 0.55, '[-1]', 'FontSize', 10, 'Interpreter', 'latex'); % Near 1 - 0.3 = 0.7
text(0.05, 1.05, '[+1]', 'FontSize', 10, 'Interpreter', 'latex'); % Near 1 + 0.3 = 1.3
text(0.05, 1.55, '[-2]', 'FontSize', 10, 'Interpreter', 'latex'); % Near 2 - 0.3 = 1.7
text(0.05, 2.05, '[+2]', 'FontSize', 10, 'Interpreter', 'latex'); % Near 2 + 0.3 = 2.3
text(0.05, 2.55, '[-3]', 'FontSize', 10, 'Interpreter', 'latex'); % Near 3 - 0.3 = 2.7
text(0.05, 3.05, '[+3]', 'FontSize', 10, 'Interpreter', 'latex'); % Near 3 + 0.3 = 3.3
text(0.05, 3.55, '[-4]', 'FontSize', 10, 'Interpreter', 'latex'); % Near 4 - 0.3 = 3.7

% Middle labels (Near the primary instability boundaries)
text(1.2, 0.35, '[-1/+0]', 'FontSize', 10, 'Interpreter', 'latex');
text(1.2, 1.35, '[-2/+1]', 'FontSize', 10, 'Interpreter', 'latex');
text(1.2, 2.35, '[-3/+2]', 'FontSize', 10, 'Interpreter', 'latex');
text(1.2, 3.35, '[-4/+3]', 'FontSize', 10, 'Interpreter', 'latex');

% --- RIGHT-SIDE LABELS (Secondary Instability Boundaries) ---
% The labels are positioned based on the value of 'w' used in the analysis.
if abs(w - 0.3) < 0.001
    % Labels for w = 0.3 
    text(4.0, 0.20, '[+0]', 'FontSize', 10, 'Interpreter', 'latex');
    text(4.0, 1.20, '[-1/+1]', 'FontSize', 10, 'Interpreter', 'latex');
    text(4.0, 2.20, '[-2/+2]', 'FontSize', 10, 'Interpreter', 'latex');
    text(4.0, 3.20, '[-3/+3]', 'FontSize', 10, 'Interpreter', 'latex');
    text(4.0, 4.20, '[-4/+4]', 'FontSize', 10, 'Interpreter', 'latex');
elseif abs(w - 0.7) < 0.001
    % Labels for w = 0.7 
    text(3.0, 0.20, '[+0]', 'FontSize', 10, 'Interpreter', 'latex');
    text(3.0, 1.20, '[-1/+1]', 'FontSize', 10, 'Interpreter', 'latex');
    text(3.0, 2.20, '[-2/+2]', 'FontSize', 10, 'Interpreter', 'latex');
    text(3.0, 3.20, '[-3/+3]', 'FontSize', 10, 'Interpreter', 'latex');
    text(3.0, 4.20, '[-4/+4]', 'FontSize', 10, 'Interpreter', 'latex');
end

% Print to Png file
print(pngfile, '-dpng')