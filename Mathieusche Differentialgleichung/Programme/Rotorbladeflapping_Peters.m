 % Damped Flapping Frequency Plot (Peters' full flapping coefficients)
%   C(t) = gamma/8 * (1 + 4*mu/3 * sin(t))
%   K(t) = p^2 + gamma/8 * (4*mu/3 * cos(t) + mu^2 * sin(2*t))
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

% --- Parameters (as in Peters) ---
p_sq = 1.0^2;          % p^2 (p = 1.0)
gamma = 12.0;          % Lock number
Omega = 1;             % fundamental (1/rev)

pngname = strrep(sprintf('Petersrotor–bladeflapping%s_%2.1f',K),'.','dot');
pngfile = fullfile(fDirPeters,[pngname,'.png']);

% base damping prefactor gamma/8
damping_pref = gamma / 8;

% --- Sweep Range (mu) ---
mu_end = 3;            % X-axis limit
mu_no = 1100;          % resolution
m_range = (-4:4);      % harmonic branches to store

mu_vals = linspace(0, mu_end, mu_no);

% prepare storage for branches
results_by_branch = struct();
for m = m_range
    if m < 0
        fname = ['m_neg_', num2str(abs(m))];
    else
        fname = ['m_', num2str(m)];
    end
    results_by_branch.(fname) = [];
end

% Floquet integration setup
T = 2*pi / Omega;
x0 = eye(2); % initial transition matrix

for k = 1:length(mu_vals)
    mu = mu_vals(k);

    % Peters' damping and stiffness as time functions:
    % C(t) = gamma/8 * (1 + (4*mu/3)*sin(t))
    % K(t) = p^2 + gamma/8 * ( (4*mu/3)*cos(t) + mu^2 * sin(2*t) )
    D_func = @(t) [0, 1;
                   -( p_sq + (gamma/8) * ( (4*mu/3)*cos(t) + (mu^2)*sin(2*t) ) ), ...
                   - (gamma/8) * (1 + (4*mu/3)*sin(t))];

    % integrate to get Phi(T) using ode45
    % reshape Phi (2x2) -> 4x1 vector for ODE
    rhs = @(t, x) reshape(D_func(t) * reshape(x, 2, 2), 4, 1);
    opts = odeset('RelTol',1e-8,'AbsTol',1e-10); % tighter tolerances help eigenvalue ordering
    [~, Phi_t] = ode45(rhs, [0, T], reshape(x0, 4, 1), opts);
    Phi_T = reshape(Phi_t(end, :), 2, 2);

    % Floquet multipliers and exponents
    Lambda = eig(Phi_T);
    eta = log(Lambda) / T;    % eta = sigma + i*omega

    normalized_omega = imag(eta) / Omega; % omega / Omega

    % branch bookkeeping: principal basis in [-0.5, 0.5], then add integers
    for r = 1:length(normalized_omega)
        freq_r = normalized_omega(r);
        basis_freq = mod(freq_r + 0.5, 1) - 0.5;
        for m = m_range
            branch_freq = basis_freq + m;
            if m < 0
                fname = ['m_neg_', num2str(abs(m))];
            else
                fname = ['m_', num2str(m)];
            end
            results_by_branch.(fname) = [results_by_branch.(fname); mu, branch_freq];
        end
    end
end

% --- Plotting ---
figure('Color','w','Position',[100 100 800 600]);
hold on;
color_map = lines(length(m_range));

title(['Frequency vs. $\mu$ (Flapping: $p=' num2str(sqrt(p_sq)) ...
       ', \gamma=' num2str(gamma) '$)'], 'Interpreter', 'latex');
xlabel('Advance Ratio ($\mu$)', 'Interpreter', 'latex');
ylabel('Frequency ($\omega/\Omega$)', 'Interpreter', 'latex');

idx = 1;
for m = m_range
    if m < 0
        fname = ['m_neg_', num2str(abs(m))];
        lbl = ['[', num2str(m), ']'];
    else
        fname = ['m_', num2str(m)];
        lbl = ['[+', num2str(m), ']'];
    end
    data = results_by_branch.(fname);
    % skip if empty (safety)
    if isempty(data)
        idx = idx + 1;
        continue;
    end

     % Plot using dots/scatters to show calculated points, which, when dense, form curves
    if useK == 1
        plot(data(:, 1), data(:, 2), '.', 'Color','k', 'MarkerSize', 8, 'DisplayName', lbl);
    else
        plot(data(:, 1), data(:, 2), '.', 'Color', color_map(idx, :), 'MarkerSize', 8, 'DisplayName', lbl);
    end
    idx = idx + 1;
end

axis([0 3 0 4.5]);
grid on;
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12, 'LineWidth', 1.0);
box on;

if ~useK
    legend('Location','northeastoutside','Interpreter','latex');
end

% annotations similar to Peters Fig.9
hold on;
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
    text(1.7, 4.2, '[-4/+4]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
    text(1.7, 3.2, '[-3/+3]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
    text(1.7, 2.2, '[-2/+2]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
    text(1.7, 1.2, '[-1/+1]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
    text(1.7, 0.2, '[+0]', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
hold off;

% Print to Png file
print(pngfile, '-dpng')