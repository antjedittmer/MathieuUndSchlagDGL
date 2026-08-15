%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Combined 3x2 figure
% Layout (3 rows, 2 columns):
%   Row 1: 3D plots using 2nd column of monodromy matrix
%          - Col 1: nu_0^2 = 0.25
%          - Col 2: nu_0^2 = 5
%   Row 2: 2D phase-plane plots (phi, phi_dot)
%          - Col 1: nu_0^2 = 0.25
%          - Col 2: nu_0^2 = 5
%   Row 3: Angle theta = unwrap(atan2(phi_dot, phi))
%          - Col 1: nu_0^2 = 0.25
%          - Col 2: nu_0^2 = 5
%
% Parameters:
%   D = 0.15
%   nu_0^2 = nu_C^2 = 0.25 and 5
%
% FIX vs. previous version: MathieuDGL's actual equations only match the
% state convention x(1)=phi, x(2)=phi_dot (same as phase_planeWindingplots.m's
% MathieuDGLsubfun, which this function is copied from) -- NOT
% x(1)=phi_dot, x(2)=phi as the old comment/Row-1 extraction assumed. The
% Row-1 extraction lines below are swapped accordingly so phi_mon/phidot_mon
% now match what MathieuDGL actually computes, and are consistent with the
% phi_pp/phidot_pp extraction already used for Rows 2-3.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
close all;

%% Global parameters
D       = 0.15;
Omega   = 1;
T       = 2*pi/Omega;
tspan   = linspace(0, T, 4000);
t0      = 0;

% Requested parameter values (nu_0^2 = nu_C^2)
nu_squared_vals = [0.25, 5];

% Initial condition for 2nd column of monodromy matrix: e2 = [0; 1]
% (state = [phi; phi_dot], so this means phi(0)=0, phi_dot(0)=1)
x0_monodromy = [0; 1];

% Initial condition for phase-plane / angle plots
x0_phase = [1; 0];

options = odeset('RelTol',1e-10,'AbsTol',1e-12);

%% Precompute solutions for both nu values
sol_struct = struct();
for idx = 1:length(nu_squared_vals)
    nu_02 = nu_squared_vals(idx);
    nu_C2 = nu_squared_vals(idx);

    % === Row 1 data: 2nd column of monodromy matrix ===
    sol_mon = ode45(@(psi,x) MathieuDGL(psi, x, D, nu_02, nu_C2), ...
                    [t0, T], x0_monodromy, options);
    y_mon   = deval(sol_mon, tspan);
    phi_mon    = y_mon(1,:)';  % phi (position)      -- FIXED (was y_mon(2,:))
    phidot_mon = y_mon(2,:)';  % phi_dot (velocity)   -- FIXED (was y_mon(1,:))

    % === Row 2 & 3 data: phase plane and angle ===
    [psi_sol, sol_raw] = ode45(@(psi,x) MathieuDGL(psi, x, D, nu_02, nu_C2), ...
                               tspan, x0_phase, options);
    phi_pp   = sol_raw(:,1); % phi
    phidot_pp = sol_raw(:,2); % phi_dot

    theta = unwrap(atan2(phidot_pp, phi_pp));

    sol_struct(idx).nu_squared   = nu_02;
    sol_struct(idx).tspan        = tspan;
    sol_struct(idx).psi_sol      = psi_sol;
    sol_struct(idx).phi_mon      = phi_mon;
    sol_struct(idx).phidot_mon   = phidot_mon;
    sol_struct(idx).phi_pp       = phi_pp;
    sol_struct(idx).phidot_pp    = phidot_pp;
    sol_struct(idx).theta        = theta;
end

%% Create combined 3x2 figure
fig = figure('Color','w', ...
    'Name','Task 1: Combined figure', ...
    'Position',[100 100 1300 900]);

for idx = 1:length(nu_squared_vals)
    nu_02 = nu_squared_vals(idx);
    col = idx;

    % --- Row 1: 3D plot (2nd column of monodromy matrix) ---
    subplot(3,2, col);
    plot3(sol_struct(idx).tspan, ...
          sol_struct(idx).phi_mon, ...
          sol_struct(idx).phidot_mon, 'b', 'LineWidth', 1.4);
    grid on;
    daspect([1 1 1]);
    xlabel('$\psi \; [rad]$','Interpreter','latex');
    ylabel('$\phi \; [-]$','Interpreter','latex');
    zlabel('$\dot{\phi} \; [-]$','Interpreter','latex');
    title(sprintf('3D plot, $\\nu_0^2=\\nu_C^2=%.2f$ (2nd monodromy col.)', nu_02), ...
          'Interpreter','latex');
    view([50 20]);

    % --- Row 2: 2D phase-plane plot (phi, phi_dot) ---
    subplot(3,2, col + 2);
    plot(sol_struct(idx).phi_pp, sol_struct(idx).phidot_pp, 'b', 'LineWidth', 1.6);
    hold on;
    plot(sol_struct(idx).phi_pp(1), sol_struct(idx).phidot_pp(1), ...
         'go','MarkerFaceColor','g','MarkerSize',8);
    plot(sol_struct(idx).phi_pp(end), sol_struct(idx).phidot_pp(end), ...
         'ro','MarkerFaceColor','r','MarkerSize',8);
    xline(0,'k--','LineWidth',1.0);
    yline(0,'k--','LineWidth',1.0);
    hold off;
    axis equal;
    grid on;
    xlabel('$\phi \; [-]$','Interpreter','latex');
    ylabel('$\dot{\phi} \; [-]$','Interpreter','latex');
    title(sprintf('Phase plane, $\\nu_0^2=\\nu_C^2=%.2f$', nu_02), ...
          'Interpreter','latex');

    % --- Row 3: Angle theta = unwrap(atan2(phi_dot, phi)) ---
    subplot(3,2, col + 4);
    plot(sol_struct(idx).psi_sol, sol_struct(idx).theta, 'm', 'LineWidth', 1.6);
    grid on;
    ax = gca;
    ax.XTick = 0:pi/2:2*pi;
    ax.XTickLabel = {'0','\pi/2','\pi','3\pi/2','2\pi'};
    xlabel('$\psi \; [rad]$','Interpreter','latex');
    ylabel('$\theta(\psi)$','Interpreter','latex');
    title(sprintf('Angle $\\theta = \\mathrm{unwrap}(\\mathrm{atan2}(\\dot\\phi,\\phi))$, $\\nu_0^2=%.2f$', nu_02), ...
          'Interpreter','latex');
end

sgtitle('Task 1: Combined winding-number visualization ($\nu_0^2 = \nu_C^2 = 0.25$ and $5$)', ...
        'Interpreter','latex','FontSize',14);

%% Save figure (optional - comment out if not needed)
% fDir = 'figureFolder';
% if ~isfolder(fDir)
%     mkdir(fDir);
% end
% fileName = fullfile(fDir, 'Task1_CombinedFigure_nu0dot25_5.svg');
% exportgraphics(fig, fileName, 'ContentType', 'vector');
% disp('Task 1: Combined 3x2 figure created and saved.');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local ODE function (Mathieu equation)
% phi'' + 2*D*phi' + (nu_0^2 + nu_C^2*cos(psi))*phi = 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dxdpsi = MathieuDGL(psi, x, D, nu_02, nu_C2)
% x(1) = phi (position)     -- FIXED comment (matches equations below)
% x(2) = phi_dot (velocity) -- FIXED comment
K_psi = nu_02 + nu_C2 * cos(psi);
dxdpsi = [x(2); -2*D*x(2) - K_psi*x(1)];
end