
% function Frequency_Peters
% Floquet Frequency Plot for Mathieu's Equation (Figure 2 Appearance)
% Separates branches according to the integer multiple 'm' used.
%
% Based on:
% David A. Peters, Sydnie M. Lieb, Loren A. Ahaus
% "Interpretation of Floquet Eigenvalues and Eigenvectors for Periodic Systems"
% JOURNAL OF THE AMERICAN HELICOPTER SOCIETY 56, 032001 (2011)
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

% --- NEW: Setup for Data Saving ---
dDir = 'dataFolder'; % Folder for Excel and .mat files
if ~isdir(dDir)
    mkdir(dDir)
end

% Plotting style selection
K = 'ColoredLines';
% K = 'BlackLines';
useK = strcmp(K,'BlackLines');
% -----------------------------------------------------------------------
% --- Parameters and Initialization ---
% -----------------------------------------------------------------------
Omega = 1;       % Fundamental frequency (Omega = 1 rad/s)
T = 2*pi/Omega; % Period
% Outer loop for different unperturbed frequencies (w)
w_values = [0.3, 0.5, 0.7];
for w = w_values
    w_sq = w^2;
    % --- Basis Frequency omega0 Calculation (Peters' Convention) ---
    basis_freq = mod(w, Omega);
    if basis_freq > Omega/2
        basis_freq = Omega - basis_freq;
    end
    omega0 = basis_freq;
    % Filename generation for saving the plot
    pngname = strrep(sprintf('PetersFrequency%s_w%1.1f',K,w),'.','dot');
    pngfile = fullfile(fDirPeters,[pngname,'.png']);
    
    % --- NEW: Filename generation for Data ---
    % Explicitly mentions "Peters" and the value of w
    dataFileName = sprintf('Peters_Mathieu_ODE_w_%1.1f', w);
    dataFileName = strrep(dataFileName, '.', 'p'); % replaces '.' with 'p' for safety
    
    % Determine the range of epsilon (x-axis limit)
    if abs(w - 0.3) < 0.001
        eps_end = 5;
        eps_no = 1000; % High resolution for w=0.3
    elseif abs(w - 0.7) < 0.001
        eps_end = 3.5; % Axis limit for w=0.7 plot
        eps_no = 150;
    else % Case for w=0.5 and others
        eps_end = 3.5;
        eps_no = 150;
    end
    eps_vals = linspace(0, eps_end, eps_no)'; % Ensure column vector for table
    m_range = (-4:4); % Integer multiple range
   
    % Structure to hold results, organized by branch 'm'
    results_by_branch = struct();
    for m = m_range
        if m < 0
            field_name = ['m_neg_', num2str(abs(m))];
        else
            field_name = ['m_', num2str(m)];
        end
        results_by_branch.(field_name) = []; 
    end
    
    % -----------------------------------------------------------------------
    % --- Floquet Exponent Calculation and Branch Separation ---
    % -----------------------------------------------------------------------
    x0 = eye(2); 
    for k = 1:length(eps_vals)
        epsilon = eps_vals(k);
        D_func = @(t) [0, 1; -(w_sq + epsilon*sin(t)), 0];
        [~, Phi_t] = ode45(@(t, x) reshape(D_func(t) * reshape(x, 2, 2), 4, 1), [0, T], reshape(x0, 4, 1));
        Phi_T = reshape(Phi_t(end, :), 2, 2); 
       
        Lambda = eig(Phi_T);
        eta = log(Lambda) / T;
        normalized_omega = imag(eta) / Omega;
        basis_freq_r = normalized_omega(1);
        
        for m = m_range
            if m==0
                branch_freq = basis_freq_r;
            else
                branch_freq = abs(m) *Omega + sign(m) *basis_freq_r;
            end
            
            if m < 0
                field_name = ['m_neg_', num2str(abs(m))];
            else
                field_name = ['m_', num2str(m)];
            end
            results_by_branch.(field_name) = [results_by_branch.(field_name); branch_freq];
        end
    end
    
    % --- NEW: Create Table for Export ---
    % The table includes Epsilon and all calculated branches
    dataTable = table(eps_vals, 'VariableNames', {'Epsilon'});
    for m = m_range
        if m < 0
            field_name = ['m_neg_', num2str(abs(m))];
        else
            field_name = ['m_', num2str(m)];
        end
        % Add each branch as a new column in the table
        dataTable.(field_name) = results_by_branch.(field_name);
    end
    
    % --- NEW: Save Excel and .mat files ---
    writetable(dataTable, fullfile(dDir, [dataFileName, '.xlsx']));
    save(fullfile(dDir, [dataFileName, '.mat']), 'dataTable');

    % -----------------------------------------------------------------------
    % --- Plotting (Original Code Continues) ---
    % -----------------------------------------------------------------------
    figure;
    hold on;
    color_map = lines;
    color_map = [color_map(1:7,:);0*ones(1,3);0.5*ones(1,3)];
    ode_str = '$\dot{x} + (w^2 + \epsilon\sin(\Omega t)) x = 0$';
    w_str = num2str(w, '%1.1f');
    new_title = {'Frequency vs. $\epsilon$, ', [ode_str, ', $w = ', w_str, '$ ($\Omega = 2\pi$ rad/s)']};
    title(new_title, 'FontSize', 16, 'Interpreter', 'latex');
    xlabel('$\epsilon$', 'FontSize', 14, 'Interpreter', 'latex');
    ylabel('Frequency ($\omega/\Omega$)', 'FontSize', 14, 'Interpreter', 'latex');
    idx = 1;
    for m = m_range
        if m < 0
            field_name = ['m_neg_', num2str(abs(m))];
        else
            field_name = ['m_', num2str(m)];
        end
        
        if m >= 0
            freq_normalized = omega0/Omega + m;
            freq_str = sprintf('$m=%+d \\rightarrow \\omega/\\Omega \\approx %.1f$', m, freq_normalized);
        elseif m < 0
            m_abs = abs(m);
            freq_normalized = m_abs - omega0/Omega;
            freq_str = sprintf('$m=%+d \\rightarrow \\omega/\\Omega \\approx %.1f$', m, freq_normalized);
        end

        y_data = results_by_branch.(field_name);
        current_color = color_map(idx, :);
        if useK == 1
            plot(eps_vals, y_data, '.', 'Color','k', 'MarkerSize', 8, 'DisplayName', freq_str);
        else
            plot(eps_vals, y_data, '.', 'Color', current_color, 'MarkerSize', 8, 'DisplayName', freq_str);
        end
        idx = idx + 1;
    end
    grid on;
    set(gca, 'TickLabelInterpreter', 'latex');
    if ~useK
        legend('Location', 'northeastoutside', 'Interpreter', 'latex');
    end
   
    % --- Branch Label Annotations (Original Code) ---
    hold on;
    text(0.05, 0.08, '[+0]', 'FontSize', 10, 'Interpreter', 'latex'); 
    text(0.05, 0.55, '[-1]', 'FontSize', 10, 'Interpreter', 'latex'); 
    text(0.05, 1.08, '[+1]', 'FontSize', 10, 'Interpreter', 'latex'); 
    text(0.05, 1.55, '[-2]', 'FontSize', 10, 'Interpreter', 'latex'); 
    text(0.05, 2.08, '[+2]', 'FontSize', 10, 'Interpreter', 'latex'); 
    text(0.05, 2.55, '[-3]', 'FontSize', 10, 'Interpreter', 'latex'); 
    text(0.05, 3.08, '[+3]', 'FontSize', 10, 'Interpreter', 'latex'); 
    text(0.05, 3.55, '[-4]', 'FontSize', 10, 'Interpreter', 'latex'); 
    text(1.2, 0.35, '[-1/+0]', 'FontSize', 10, 'Interpreter', 'latex');
    text(1.2, 1.35, '[-2/+1]', 'FontSize', 10, 'Interpreter', 'latex');
    text(1.2, 2.35, '[-3/+2]', 'FontSize', 10, 'Interpreter', 'latex');
    text(1.2, 3.35, '[-4/+3]', 'FontSize', 10, 'Interpreter', 'latex');
    if abs(w - 0.3) < 0.001
        text(4.0, 0.20, '[+0]', 'FontSize', 10, 'Interpreter', 'latex');
        text(4.0, 1.20, '[-1/+1]', 'FontSize', 10, 'Interpreter', 'latex');
        text(4.0, 2.20, '[-2/+2]', 'FontSize', 10, 'Interpreter', 'latex');
        text(4.0, 3.20, '[-3/+3]', 'FontSize', 10, 'Interpreter', 'latex');
        text(4.0, 4.20, '[-4/+4]', 'FontSize', 10, 'Interpreter', 'latex');
    elseif abs(w - 0.7) < 0.001
        text(3.0, 0.20, '[+0]', 'FontSize', 10, 'Interpreter', 'latex');
        text(3.0, 1.20, '[-1/+1]', 'FontSize', 10, 'Interpreter', 'latex');
        text(3.0, 2.20, '[-2/+2]', 'FontSize', 10, 'Interpreter', 'latex');
        text(3.0, 3.20, '[-3/+3]', 'FontSize', 10, 'Interpreter', 'latex');
        text(3.0, 4.20, '[-4/+4]', 'FontSize', 10, 'Interpreter', 'latex');
    end
    print(pngfile, '-dpng')
=======
% function Frequency_Peters
% Floquet Frequency Plot for Mathieu's Equation (Figure 2 Appearance)
% Separates branches according to the integer multiple 'm' used.
%
% Based on:
% David A. Peters, Sydnie M. Lieb, Loren A. Ahaus
% "Interpretation of Floquet Eigenvalues and Eigenvectors for Periodic Systems"
% JOURNAL OF THE AMERICAN HELICOPTER SOCIETY 56, 032001 (2011)
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

% --- NEW: Setup for Data Saving ---
dDir = 'dataFolder'; % Folder for Excel and .mat files
if ~isdir(dDir)
    mkdir(dDir)
end

% Plotting style selection
K = 'ColoredLines';
% K = 'BlackLines';
useK = strcmp(K,'BlackLines');
% -----------------------------------------------------------------------
% --- Parameters and Initialization ---
% -----------------------------------------------------------------------
Omega = 1;       % Fundamental frequency (Omega = 1 rad/s)
T = 2*pi/Omega; % Period
% Outer loop for different unperturbed frequencies (w)
w_values = [0.3, 0.5, 0.7];
for w = w_values
    w_sq = w^2;
    % --- Basis Frequency omega0 Calculation (Peters' Convention) ---
    basis_freq = mod(w, Omega);
    if basis_freq > Omega/2
        basis_freq = Omega - basis_freq;
    end
    omega0 = basis_freq;
    % Filename generation for saving the plot
    pngname = strrep(sprintf('PetersFrequency%s_w%1.1f',K,w),'.','dot');
    pngfile = fullfile(fDirPeters,[pngname,'.png']);
    
    % --- NEW: Filename generation for Data ---
    % Explicitly mentions "Peters" and the value of w
    dataFileName = sprintf('Peters_Mathieu_ODE_w_%1.1f', w);
    dataFileName = strrep(dataFileName, '.', 'p'); % replaces '.' with 'p' for safety
    
    % Determine the range of epsilon (x-axis limit)
    if abs(w - 0.3) < 0.001
        eps_end = 5;
        eps_no = 1000; % High resolution for w=0.3
    elseif abs(w - 0.7) < 0.001
        eps_end = 3.5; % Axis limit for w=0.7 plot
        eps_no = 150;
    else % Case for w=0.5 and others
        eps_end = 3.5;
        eps_no = 150;
    end
    eps_vals = linspace(0, eps_end, eps_no)'; % Ensure column vector for table
    m_range = (-4:4); % Integer multiple range
   
    % Structure to hold results, organized by branch 'm'
    results_by_branch = struct();
    for m = m_range
        if m < 0
            field_name = ['m_neg_', num2str(abs(m))];
        else
            field_name = ['m_', num2str(m)];
        end
        results_by_branch.(field_name) = []; 
    end
    
    % -----------------------------------------------------------------------
    % --- Floquet Exponent Calculation and Branch Separation ---
    % -----------------------------------------------------------------------
    x0 = eye(2); 
    for k = 1:length(eps_vals)
        epsilon = eps_vals(k);
        D_func = @(t) [0, 1; -(w_sq + epsilon*sin(t)), 0];
        [~, Phi_t] = ode45(@(t, x) reshape(D_func(t) * reshape(x, 2, 2), 4, 1), [0, T], reshape(x0, 4, 1));
        Phi_T = reshape(Phi_t(end, :), 2, 2); 
       
        Lambda = eig(Phi_T);
        eta = log(Lambda) / T;
        normalized_omega = imag(eta) / Omega;
        basis_freq_r = normalized_omega(1);
        
        for m = m_range
            if m==0
                branch_freq = basis_freq_r;
            else
                branch_freq = abs(m) *Omega + sign(m) *basis_freq_r;
            end
            
            if m < 0
                field_name = ['m_neg_', num2str(abs(m))];
            else
                field_name = ['m_', num2str(m)];
            end
            results_by_branch.(field_name) = [results_by_branch.(field_name); branch_freq];
        end
    end
    
    % --- NEW: Create Table for Export ---
    % The table includes Epsilon and all calculated branches
    dataTable = table(eps_vals, 'VariableNames', {'Epsilon'});
    for m = m_range
        if m < 0
            field_name = ['m_neg_', num2str(abs(m))];
        else
            field_name = ['m_', num2str(m)];
        end
        % Add each branch as a new column in the table
        dataTable.(field_name) = results_by_branch.(field_name);
    end
    
    % --- NEW: Save Excel and .mat files ---
    writetable(dataTable, fullfile(dDir, [dataFileName, '.xlsx']));
    save(fullfile(dDir, [dataFileName, '.mat']), 'dataTable');

    % -----------------------------------------------------------------------
    % --- Plotting (Original Code Continues) ---
    % -----------------------------------------------------------------------
    figure;
    hold on;
    color_map = lines;
    color_map = [color_map(1:7,:);0*ones(1,3);0.5*ones(1,3)];
    ode_str = '$\dot{x} + (w^2 + \epsilon\sin(\Omega t)) x = 0$';
    w_str = num2str(w, '%1.1f');
    new_title = {'Frequency vs. $\epsilon$, ', [ode_str, ', $w = ', w_str, '$ ($\Omega = 2\pi$ rad/s)']};
    title(new_title, 'FontSize', 16, 'Interpreter', 'latex');
    xlabel('$\epsilon$', 'FontSize', 14, 'Interpreter', 'latex');
    ylabel('Frequency ($\omega/\Omega$)', 'FontSize', 14, 'Interpreter', 'latex');
    idx = 1;
    for m = m_range
        if m < 0
            field_name = ['m_neg_', num2str(abs(m))];
        else
            field_name = ['m_', num2str(m)];
        end
        
        if m >= 0
            freq_normalized = omega0/Omega + m;
            freq_str = sprintf('$m=%+d \\rightarrow \\omega/\\Omega \\approx %.1f$', m, freq_normalized);
        elseif m < 0
            m_abs = abs(m);
            freq_normalized = m_abs - omega0/Omega;
            freq_str = sprintf('$m=%+d \\rightarrow \\omega/\\Omega \\approx %.1f$', m, freq_normalized);
        end

        y_data = results_by_branch.(field_name);
        current_color = color_map(idx, :);
        if useK == 1
            plot(eps_vals, y_data, '.', 'Color','k', 'MarkerSize', 8, 'DisplayName', freq_str);
        else
            plot(eps_vals, y_data, '.', 'Color', current_color, 'MarkerSize', 8, 'DisplayName', freq_str);
        end
        idx = idx + 1;
    end
    grid on;
    set(gca, 'TickLabelInterpreter', 'latex');
    if ~useK
        legend('Location', 'northeastoutside', 'Interpreter', 'latex');
    end
   
    % --- Branch Label Annotations (Original Code) ---
    hold on;
    text(0.05, 0.08, '[+0]', 'FontSize', 10, 'Interpreter', 'latex'); 
    text(0.05, 0.55, '[-1]', 'FontSize', 10, 'Interpreter', 'latex'); 
    text(0.05, 1.08, '[+1]', 'FontSize', 10, 'Interpreter', 'latex'); 
    text(0.05, 1.55, '[-2]', 'FontSize', 10, 'Interpreter', 'latex'); 
    text(0.05, 2.08, '[+2]', 'FontSize', 10, 'Interpreter', 'latex'); 
    text(0.05, 2.55, '[-3]', 'FontSize', 10, 'Interpreter', 'latex'); 
    text(0.05, 3.08, '[+3]', 'FontSize', 10, 'Interpreter', 'latex'); 
    text(0.05, 3.55, '[-4]', 'FontSize', 10, 'Interpreter', 'latex'); 
    text(1.2, 0.35, '[-1/+0]', 'FontSize', 10, 'Interpreter', 'latex');
    text(1.2, 1.35, '[-2/+1]', 'FontSize', 10, 'Interpreter', 'latex');
    text(1.2, 2.35, '[-3/+2]', 'FontSize', 10, 'Interpreter', 'latex');
    text(1.2, 3.35, '[-4/+3]', 'FontSize', 10, 'Interpreter', 'latex');
    if abs(w - 0.3) < 0.001
        text(4.0, 0.20, '[+0]', 'FontSize', 10, 'Interpreter', 'latex');
        text(4.0, 1.20, '[-1/+1]', 'FontSize', 10, 'Interpreter', 'latex');
        text(4.0, 2.20, '[-2/+2]', 'FontSize', 10, 'Interpreter', 'latex');
        text(4.0, 3.20, '[-3/+3]', 'FontSize', 10, 'Interpreter', 'latex');
        text(4.0, 4.20, '[-4/+4]', 'FontSize', 10, 'Interpreter', 'latex');
    elseif abs(w - 0.7) < 0.001
        text(3.0, 0.20, '[+0]', 'FontSize', 10, 'Interpreter', 'latex');
        text(3.0, 1.20, '[-1/+1]', 'FontSize', 10, 'Interpreter', 'latex');
        text(3.0, 2.20, '[-2/+2]', 'FontSize', 10, 'Interpreter', 'latex');
        text(3.0, 3.20, '[-3/+3]', 'FontSize', 10, 'Interpreter', 'latex');
        text(3.0, 4.20, '[-4/+4]', 'FontSize', 10, 'Interpreter', 'latex');
    end
    print(pngfile, '-dpng')
>>>>>>> 16342c24dd3d857b752736907883e637769b53ac
end