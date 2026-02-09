% Damped Flapping Frequency Plot (Peters' full flapping coefficients)
%   C(t) = gamma/8 * (1 + 4*mu/3 * sin(t))
%   K(t) = p^2 + gamma/8 * (4*mu/3 * cos(t) + mu^2 * sin(2*t))
clear; clc; close all;
fDir = 'figureFolder'; % Folder for figures
if ~isdir(fDir)
    mkdir(fDir)
end
fDirPeters = fullfile(fDir,'figureFolderPeters'); % Subfolder specific to Peters' plots
if ~isdir(fDirPeters)
    mkdir(fDirPeters)
end
K = 'ColoredLines';
useK = strcmp(K,'BlackLines');

% --- Parameters ---
p_sq = 1.0^2;          
gamma_list = [12.0, 9.6];  
Omega = 1;             
mu_end = 3;            
mu_no = 1100;              
mu_vals = linspace(0, mu_end, mu_no);
T = 2*pi / Omega;
x0 = eye(2); 

% --- Outer Loop for different Gamma values ---
for g_val = gamma_list
    gamma = g_val;
    fprintf('Processing gamma = %.1f...\n', gamma);

    % *** GAMMA-SPECIFIC HARMONIC RANGE ***
    if gamma == 12.0
        m_range = (-4:4);  % 9 branches for gamma=12
    else  % gamma=9.6
        m_range = (-3:3);  % 7 branches for gamma=9.6
    end
    
    % prepare storage for branches for this specific gamma
    results_by_branch = struct();
    for m = m_range
        fname = matlab.lang.makeValidName(['m_' num2str(m)]);
        results_by_branch.(fname) = [];
    end
    
    % Floquet integration sweep
    for k = 1:length(mu_vals)
        mu = mu_vals(k);
        
        % Peters' system matrix
        D_func = @(t) [0, 1;
                       -( p_sq + (gamma/8) * ( (4*mu/3)*cos(t) + (mu^2)*sin(2*t) ) ), ...
                       - (gamma/8) * (1 + (4*mu/3)*sin(t))];
        
        rhs = @(t, x) reshape(D_func(t) * reshape(x, 2, 2), 4, 1);
        opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
        
        [~, Phi_t] = ode45(rhs, [0, T], reshape(x0, 4, 1), opts);
        Phi_T = reshape(Phi_t(end, :), 2, 2);
        
        Lambda = eig(Phi_T);
        eta = log(Lambda) / T;    
        normalized_omega = imag(eta) / Omega; 
        
        for r = 1:length(normalized_omega)
            freq_r = normalized_omega(r);
            basis_freq = mod(freq_r + 0.5, 1) - 0.5;
            for m = m_range
                branch_freq = basis_freq + m;
                fname = matlab.lang.makeValidName(['m_' num2str(m)]);
                results_by_branch.(fname) = [results_by_branch.(fname); mu, branch_freq];
            end
        end
    end
    
    % --- Plotting for this Gamma ---
    fig = figure('Color','w','Name', sprintf('Gamma %.1f - Peters Flapping', gamma), ...
                 'Position',[100 100 1200 600]);
    hold on;
    
    % Create consistent color map for all 9 harmonics (-4:4)
    colors = lines(length(m_range));  
    
    title(['Frequency vs. $\mu$ (Flapping: $p=' num2str(sqrt(p_sq)) ...
           ', \gamma=' num2str(gamma) '$)'], 'Interpreter', 'latex');
    xlabel('$\mu$', 'Interpreter', 'latex');
    ylabel('Frequency ($\omega/\Omega$)', 'Interpreter', 'latex');
    
    % Single clean loop - each m gets unique color by index
    for i = 1:length(m_range)
        m = m_range(i);
        fname = matlab.lang.makeValidName(['m_' num2str(m)]);
        data = results_by_branch.(fname);
        
        if isempty(data), continue; end
        
        plot_color = colors(i,:);
        if useK, plot_color = 'k'; end  
        
        plot(data(:,1), data(:,2), '.', 'Color', plot_color, ...
             'MarkerSize', 8, 'DisplayName', sprintf('$m=%+d$', m));
    end
    
    % Gamma-specific harmonic labels with updated axis limits
    if gamma == 12.0
        % Labels for gamma = 12.0 (top-right corner at x=3, y=4.5)
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
        axis([0 3 0 4.5]);  
    else  % gamma = 9.6
        % Specific labels for gamma = 9.6
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

        axis([0 2 0 3.5]);  
    end
    
    grid on; box on;
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
    
    if ~useK
        legend('Location','northeastoutside','Interpreter','latex');
    end
    
    % Save result
    pngname = strrep(sprintf('Peters_RotorbladeFlapping_%2.1f', gamma), '.', 'dot');
    print(fullfile(fDirPeters, [pngname '.png']), '-dpng');
end
