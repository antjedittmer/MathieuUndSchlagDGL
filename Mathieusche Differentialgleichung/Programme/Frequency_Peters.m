% MATLAB Code to Recreate Mathieu's Equation Frequency Plot (Figure 2 Appearance)
% Separates branches according to the integer multiple 'm' used.

clear; clc; close all;

fDir = 'figureFolder'; % Ordner Abbildungen
if ~isdir(fDir) %#ok<ISDIR>
    mkdir(fDir)
end

fDirPeters = fullfile(fDir,'figureFolderPeters');
if ~isdir(fDirPeters) %#ok<ISDIR>
    mkdir(fDirPeters)
end

K = 'ColoredLines';
% K = 'BlackLines';
useK = strcmp(K,'BlackLines');

% --- Parameters ---
w_sq = 0.7^2;    % Natural frequency squared (w=0.7, as in the paper)
w = sqrt(w_sq);  % Natural frequency w = 0.7
Omega = 1;       % Fundamental frequency (1/rev)

pngname = strrep(sprintf('PetersFrequency%s_%2.1f.png',K,w),'.','dot');
pngfile = fullfile(fDirPeters,pngname);


% Updated range to cover requested x-axis limit (0 to 3.5)
if abs(w - 0.3) <0.001
    eps_end = 5; %  x axis limit in AHS2009 paper
    eps_no = 1000; % needs to be high to visualize steep flank at 3.5
else
    eps_end = 3.5; % axis limit in JAHS2011 paper
    eps_no = 150; % original value
end

eps_vals = linspace(0, eps_end, eps_no);
m_range = (-4:4); % Integer multiple range for plotting branches

% Structure to hold results, organized by branch 'm'
results_by_branch = struct();
for m = m_range
    % FIX: Create a valid field name (e.g., 'm_neg_4', 'm_0', 'm_1')
    if m < 0
        field_name = ['m_neg_', num2str(abs(m))];
    else
        field_name = ['m_', num2str(m)];
    end
    results_by_branch.(field_name) = [];
end

% --- Floquet Exponent Calculation ---
T = 2*pi/Omega; % Period
x0 = eye(2); % Intial state as identity matrix

for k = 1:length(eps_vals)
    epsilon = eps_vals(k);

    % The state matrix D(t) for x' = D(t)x, where x = [x; x_dot]
    D_func = @(t) [0, 1; -(w_sq + epsilon*sin(t)), 0];


    % Solve the Transition Matrix Phi(T)
    [~, Phi_t] = ode45(@(t, x) reshape(D_func(t) * reshape(x, 2, 2), 4, 1), [0, T], reshape(x0, 4, 1));
    Phi_T = reshape(Phi_t(end, :), 2, 2); % Transition Matrix at T

    % Floquet Exponents (eta = 1/T * log(Lambda))
    Lambda = eig(Phi_T);
    eta = log(Lambda) / T;

    % Extract Imaginary Part (normalized frequency), omega/Omega
    normalized_omega = imag(eta) / Omega;

    % --- Separate and Store Branches (m) ---
    for r = 1:length(normalized_omega)
        freq_r = normalized_omega(r);

        % 1. Find the basis frequency (principal value in the range [-0.5, 0.5])
        basis_freq = mod(freq_r + 0.5, 1) - 0.5;

        % 2. Calculate and store all possible branches m for this single root 'r'
        for m = m_range
            branch_freq = basis_freq + m;

            % FIX: Determine the valid field name
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

% --- Plotting ---
figure;
hold on;
color_map = lines(length(m_range)); % Use a color map for distinct lines

% Title, XLabel, and YLabel
title(['Frequency vs. $\epsilon$ (for $w = ' num2str(w) '$)'], 'Interpreter', 'latex');
xlabel('$\epsilon$', 'Interpreter', 'latex');
ylabel('Frequency ($\omega$)', 'Interpreter', 'latex');

m_labels = {};
idx = 1;
for m = m_range
    % Determine field name and label string
    if m < 0
        field_name = ['m_neg_', num2str(abs(m))];
        label = ['[', num2str(m), ']']; % e.g., [-1]
    else
        field_name = ['m_', num2str(m)];
        label = ['[+', num2str(m), ']']; % e.g., [+0]
    end

    % Plot the data for this specific branch 'm'
    data = results_by_branch.(field_name);

    % Plot using dots/scatters to show calculated points, which, when dense, form curves
    if useK == 1
        plot(data(:, 1), data(:, 2), '.', 'Color','k', 'MarkerSize', 8, 'DisplayName', label);
    else
        plot(data(:, 1), data(:, 2), '.', 'Color', color_map(idx, :), 'MarkerSize', 8, 'DisplayName', label);
    end
    idx = idx + 1;
end

% Set Axis limits as requested
axis([0 eps_end 0 4.5]);
grid on;

set(gca, 'TickLabelInterpreter', 'latex');
hold off;
% --- Add branch labels exactly as in the figure ---
hold on;

% Left-side labels
text(0.05, 0.05, '[+0]', 'FontSize', 10, 'Interpreter', 'latex');
text(0.05, 0.55, '[-1]', 'FontSize', 10, 'Interpreter', 'latex');
text(0.05, 1.05, '[+1]', 'FontSize', 10, 'Interpreter', 'latex');
text(0.05, 1.55, '[-2]', 'FontSize', 10, 'Interpreter', 'latex');
text(0.05, 2.05, '[+2]', 'FontSize', 10, 'Interpreter', 'latex');
text(0.05, 2.55, '[-3]', 'FontSize', 10, 'Interpreter', 'latex');
text(0.05, 3.05, '[+3]', 'FontSize', 10, 'Interpreter', 'latex');
text(0.05, 3.55, '[-4]', 'FontSize', 10, 'Interpreter', 'latex');

% Middle labels (ε ≈ 1–2)
text(1.2, 0.35, '[-1/+0]', 'FontSize', 10, 'Interpreter', 'latex');
text(1.2, 1.35, '[-2/+1]', 'FontSize', 10, 'Interpreter', 'latex');
text(1.2, 2.35, '[-3/+2]', 'FontSize', 10, 'Interpreter', 'latex');
text(1.2, 3.35, '[-4/+3]', 'FontSize', 10, 'Interpreter', 'latex');

% Right-side labels (ε ≈ 2.7–3)
text(2.8, 0.05, '[+0]', 'FontSize', 10, 'Interpreter', 'latex');
text(2.8, 0.55, '[-1/+1]', 'FontSize', 10, 'Interpreter', 'latex');
text(2.8, 1.55, '[-2/+2]', 'FontSize', 10, 'Interpreter', 'latex');
text(2.8, 2.55, '[-3/+3]', 'FontSize', 10, 'Interpreter', 'latex');
text(2.8, 3.55, '[-4/+4]', 'FontSize', 10, 'Interpreter', 'latex');

hold off;

% Print to Png file
print(pngfile, '-dpng')