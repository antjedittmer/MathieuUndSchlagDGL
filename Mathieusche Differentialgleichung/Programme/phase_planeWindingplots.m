%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D winding-number visualization in phase plane
%
% Based on the winding-number logic from:
% MathieuStruttscheKarte_Peters.m
%
% Equation:
%   phi'' + 2D*phi' + (nu_0^2 + nu_C^2*cos(psi))*phi = 0
%
% For the requested cases:
%   nu_0^2 = nu_C^2 = 0.25, 5, 8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear;
close all;

%% Parameters
D = 0.15;              % consistent with Peters script
Omega = 1;
T = 2*pi/Omega;
tspan = linspace(0, T, 4000);

% Requested values
nu_vals = [0.25, 5, 8];

% Initial condition for visualization
x0 = [1; 0];

%% Figure folders
fDir = 'figureFolder';
if ~isfolder(fDir)
    mkdir(fDir);
end

fDirTask2 = fullfile(fDir, 'Winding2D_SVG');
if ~isfolder(fDirTask2)
    mkdir(fDirTask2);
end

%% Loop over requested nu values
for idx = 1:length(nu_vals)

    nu_02 = nu_vals(idx);
    nu_C2 = nu_vals(idx);

    options = odeset('RelTol',1e-10,'AbsTol',1e-12);

    % Solve system
    % x(1) = phi
    % x(2) = phi_dot
    [psi_sol, sol_raw] = ode45(@(psi, x) MathieuDGLsubfun(psi, x, D, nu_02, nu_C2), ...
                               tspan, x0, options);

    x_path     = sol_raw(:,1);   % phi
    x_dot_path = sol_raw(:,2);   % phi_dot

    % Winding number calculation
    theta = unwrap(atan2(x_dot_path, x_path));
    total_delta_theta = abs(theta(end) - theta(1));
    m_direct = round(total_delta_theta / pi) * 0.5;

    %% ==============================================================
    % Figure 1: phase plane (phi, phi_dot)
    % ==============================================================
    fig1 = figure('Color','w', ...
                  'Name', sprintf('Phase plane nu = %.2f', nu_02), ...
                  'Position', [100 100 850 700]);
    hold on;
    grid on;
    box on;

    plot(x_path, x_dot_path, 'b', 'LineWidth', 1.8);
    plot(x_path(1),   x_dot_path(1),   'go', 'MarkerFaceColor','g', 'MarkerSize',8);
    plot(x_path(end), x_dot_path(end), 'ro', 'MarkerFaceColor','r', 'MarkerSize',8);

    xline(0, 'k--', 'LineWidth', 1.0);
    yline(0, 'k--', 'LineWidth', 1.0);

    xlabel('$\phi \; [-]$', 'Interpreter','latex');
    ylabel('$\dot{\phi} \; [-]$', 'Interpreter','latex');
    title(sprintf('Phase-plane trajectory for $\\nu_0^2 = \\nu_C^2 = %.2f$, $D = %.2f$', nu_02, D), ...
          'Interpreter','latex');

    axis equal;
    ax = gca;
    ax.FontSize = 12;

    % Safe text annotation
    xl = xlim;
    yl = ylim;
    txt1 = sprintf('m_{direct} = %.1f', m_direct);
    text(xl(1) + 0.05*(xl(2)-xl(1)), yl(2) - 0.08*(yl(2)-yl(1)), txt1, ...
        'Interpreter','tex', 'FontSize', 13, 'BackgroundColor','w', ...
        'EdgeColor',[0.8 0.8 0.8], 'Margin', 4);

    legend({'trajectory','start','end','\phi = 0','\phi-dot = 0'}, ...
           'Interpreter','tex', 'Location','best');

    baseName1 = sprintf('Task2_PhasePlane_D%0.2f_nu%0.2f', D, nu_02);
    baseName1 = strrep(baseName1, '.', 'dot');
    fileName1 = fullfile(fDirTask2, [baseName1 '.svg']);
    exportgraphics(fig1, fileName1, 'Resolution', 300);

    %% ==============================================================
    % Figure 2: unwrapped angle theta(psi)
    % ==============================================================
    fig2 = figure('Color','w', ...
                  'Name', sprintf('Unwrapped angle nu = %.2f', nu_02), ...
                  'Position', [150 120 850 500]);
    hold on;
    grid on;
    box on;

    plot(psi_sol, theta, 'm', 'LineWidth', 1.6);

    xlabel('$\psi \; [rad]$', 'Interpreter','latex');
    ylabel('$\theta(\psi)$', 'Interpreter','latex');
    title(sprintf('Unwrapped phase angle for $\\nu_0^2 = \\nu_C^2 = %.2f$', nu_02), ...
          'Interpreter','latex');

    ax = gca;
    ax.XTick = 0:pi/2:2*pi;
    ax.XTickLabel = {'0','\pi/2','\pi','3\pi/2','2\pi'};
    ax.FontSize = 12;

    % Corrected annotation syntax
    txt2 = sprintf('\\Delta\\theta = %.3f,  m = %.1f', total_delta_theta, m_direct);
    text(psi_sol(round(end*0.08)), theta(round(end*0.90)), txt2, ...
        'Interpreter','tex', 'FontSize', 12, 'BackgroundColor','w', ...
        'EdgeColor',[0.8 0.8 0.8], 'Margin', 4);

    baseName2 = sprintf('Task2_Angle_D%0.2f_nu%0.2f', D, nu_02);
    baseName2 = strrep(baseName2, '.', 'dot');
    fileName2 = fullfile(fDirTask2, [baseName2 '.svg']);
    exportgraphics(fig2, fileName2, 'Resolution', 300);

end

disp('2D phase-plane winding plots created.');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local ODE function for task 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dxdpsi = MathieuDGLsubfun(psi, x, D, nu_02, nu_C2)
% x(1) = phi
% x(2) = phi_dot
K_psi = nu_02 + nu_C2 * cos(psi);
dxdpsi = [x(2); -2*D*x(2) - K_psi*x(1)];
end