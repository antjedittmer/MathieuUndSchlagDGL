% Damped Flapping Frequency Plot (Peters' full flapping coefficients)
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

% --- Parameters ---
p_sq = 1.0^2;
gamma_list = [12.0, 9.6];
Omega = 1;
mu_no = 1100;
T = 2*pi / Omega;
x0 = eye(2);
for g_val = gamma_list
    gamma = g_val;
    fprintf('Processing gamma = %.1f...\n', gamma);
    % Define range based on paper figures
    if gamma == 12.0
        mu_end = 3; m_range = (-4:4);
    else
        mu_end = 2; m_range = (-3:3);
    end
    mu_vals = linspace(0, mu_end, mu_no);

    results_by_branch = struct();
    for m = m_range
        fname = matlab.lang.makeValidName(['m_' num2str(m)]);
        results_by_branch.(fname) = zeros(length(mu_vals), 2);
    end

    for k = 1:length(mu_vals)
        mu = mu_vals(k);

        % Dynamic Matrix
        D_func = @(t) [0, 1;
            -( p_sq + (gamma/8) * ( (4*mu/3)*cos(t) + (mu^2)*sin(2*t) ) ), ...
            - (gamma/8) * (1 + (4*mu/3)*sin(t))];

        [~, Phi_t] = ode45(@(t, x) reshape(D_func(t) * reshape(x, 2, 2), 4, 1), ...
            [0, T], reshape(x0, 4, 1), odeset('RelTol',1e-8));
        Phi_T = reshape(Phi_t(end, :), 2, 2);

        Lambda = eig(Phi_T);
        eta = log(Lambda) / T;

        % Use absolute imaginary part of one eigenvalue
        % This prevents double-plotting and ensures the basis is positive
        basis_freq = abs(imag(eta(1))) / Omega;

        for m = m_range
            % Proper Peters branch interpretation
            % Frequencies are shifted by m and follow the sign of the harmonic
            branch_freq = abs(m + basis_freq);
            if m < 0, branch_freq = abs(abs(m) - basis_freq); end

            fname = matlab.lang.makeValidName(['m_' num2str(m)]);
            results_by_branch.(fname)(k, :) = [mu, branch_freq];
        end
    end

    % --- Plotting ---
    aFig = figure('Color','w','Name', sprintf('Gamma %.1f', gamma));
    pos0 = get(0,'defaultFigurePosition'); 
    aFig.Position = [pos0(1:2), 1.4*pos0(3), pos0(4)];
   
    hold on;
    % Define the color map: first 7 from 'lines', then black, then gray
    base_colors = lines(7);
    color_map = [base_colors; 0 0 0; 0.5 0.5 0.5];

    % If m_range is larger than our manual map, expand it to avoid indexing errors
    if length(m_range) > size(color_map, 1)
        color_map = lines(length(m_range));
    end

    for i = 1:length(m_range)
        m = m_range(i);
        fname = matlab.lang.makeValidName(['m_' num2str(m)]);
        data = results_by_branch.(fname);
        % Use the corrected variable name: color_map
        plot(data(:,1), data(:,2), '.', 'Color', color_map(i,:), 'MarkerSize', 6, ...
            'DisplayName', sprintf('m = %+d', m));
    end

    % Finalizing plot
    title(['Frequency branches, rotor blade flapping, $p=1.0, \gamma=' num2str(gamma) '$'], ...
        'Interpreter', 'latex', 'FontSize', 14);
    xlabel('$\mu$', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('Frequency $\omega/\Omega$', 'Interpreter', 'latex', 'FontSize', 12);
    grid on; box on;

    % Add Legend
    legend('Location', 'northeastoutside', 'Interpreter', 'latex');

    if gamma == 12.0, axis([0 3 0 4.5]); else, axis([0 2 0 3.5]); end
    % Add gamma-specific textual labels
    hold on;
    if abs(gamma - 12.0) < 1e-6
        % gamma = 12.0 labels
        text(0.07, 3.8, '[-4]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
        text(0.07, 3.3,  '[+3]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
        text(0.54, 3.8,'[-4/+3]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
        text(0.07, 2.8, '[-3]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
        text(0.07, 2.3,  '[+2]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
        text(0.54, 2.8, '[-3/+2]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
        text(0.07, 1.8, '[-2]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
        text(0.07, 1.3,  '[+1]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
        text(0.54, 1.8, '[-2/+1]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
        text(0.07, 0.8, '[-1]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
        text(0.07, 0.3,  '[+0]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
        text(0.54, 0.8, '[-1/+0]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
        text(2.7, 4.3, '[-4/+4]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
        text(2.7, 3.3, '[-3/+3]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
        text(2.7, 2.3, '[-2/+2]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
        text(2.7, 1.3, '[-1/+1]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
        text(2.7, 0.3, '[+0]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
    elseif abs(gamma - 9.6) < 1e-6
        % gamma = 9.6 labels
        text(1.4, 3.2, '$[-3/+3]$', 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
        text(1.4, 2.2, '$[-2/+2]$', 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
        text(1.4, 1.2, '$[-1/+1]$', 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
        text(0.4, 0.2, '$[+0]$', 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
        text(0.6, 0.80, '$[-1]$', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
        text(0.4, 1.22, '$[+1]$', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
        text(0.6, 1.80, '$[-2]$', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
        text(0.4, 2.22, '$[+2]$', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
        text(0.6, 2.80, '$[-3]$', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
        text(0.4, 3.22, '$[+3]$', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
    end
    
    % --- File Naming and Saving ---
    gStr = strrep(num2str(gamma), '.', 'p');
    pngname = sprintf('PetersRotorFlapping_p1p0_gamma%s', gStr);
    pngfile = fullfile(fDirPeters,[pngname,'.png']);
    print(pngfile, '-dpng')
end