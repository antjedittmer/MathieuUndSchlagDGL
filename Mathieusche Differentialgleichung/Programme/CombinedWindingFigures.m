%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Combined 3x2 figure
%
% Row 1: 3D plots corresponding to the second transition-matrix column
% Row 2: Phase-plane plots
% Row 3: Unwrapped phase angle
%
% nu_0^2 = nu_C^2 = 0.25 and 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear;
close all;

%% Global parameters
D       = 0.15;
Omega   = 1;
T       = 2*pi/Omega;
t0      = 0;

% Use the same discretization as MathieuTransitionMatrix.m
tspan = t0:0.0001:T;

% Parameter values
nu_squared_vals = [0.25, 5];

% Initial conditions
% First column: e1 = [1; 0]
% Second column: e2 = [0; 1]
e1 = [1; 0];
e2 = [0; 1];

% Initial condition for the phase-plane and angle plots
x0_phase = [1; 0];

options = odeset( ...
    'RelTol', 1e-10, ...
    'AbsTol', 1e-12);

%% Preallocate solution structure
sol_struct = struct();

%% Solve the equations
for idx = 1:length(nu_squared_vals)

    nu_02 = nu_squared_vals(idx);
    nu_C2 = nu_02;

    %% First transition-matrix column
    sol1 = ode45( ...
        @(psi,x) MathieuDGL(psi, x, D, nu_02, nu_C2), ...
        [t0, T], e1, options);

    y1 = deval(sol1, tspan);

    %% Second transition-matrix column
    sol2 = ode45( ...
        @(psi,x) MathieuDGL(psi, x, D, nu_02, nu_C2), ...
        [t0, T], e2, options);

    y2 = deval(sol2, tspan);

    % Transition matrix at each psi:
    % Phi(psi) = [y1(:,psi), y2(:,psi)]
    %
    % For the required first-row 3D plots, use the second column:
    % y2(1,:) = phi
    % y2(2,:) = phi_dot

    phi_mon    = y2(1, :)';   % phi, second transition-matrix column
    phidot_mon = y2(2, :)';   % phi_dot, second transition-matrix column

    %% Phase-plane solution
    [psi_sol, sol_raw] = ode45( ...
        @(psi,x) MathieuDGL(psi, x, D, nu_02, nu_C2), ...
        tspan, x0_phase, options);

    phi_pp    = sol_raw(:,1);
    phidot_pp = sol_raw(:,2);

    %% Angle
    theta = unwrap(atan2(phidot_pp, phi_pp));

    %% Monodromy matrix at psi = T
    Monodromy = [y1(:,end), y2(:,end)];

    %% Store results
    sol_struct(idx).nu_squared = nu_02;
    sol_struct(idx).tspan      = tspan;

    sol_struct(idx).phi_mon    = phi_mon;
    sol_struct(idx).phidot_mon = phidot_mon;

    sol_struct(idx).phi_pp     = phi_pp;
    sol_struct(idx).phidot_pp  = phidot_pp;
    sol_struct(idx).psi_sol    = psi_sol;
    sol_struct(idx).theta      = theta;

    sol_struct(idx).Monodromy  = Monodromy;
end

%% Create combined 3x2 figure
fig = figure( ...
    'Color', 'w', ...
    'Name', 'Combined winding-number figure', ...
    'Position', [100 100 1300 900]);

for idx = 1:length(nu_squared_vals)

    nu_02 = sol_struct(idx).nu_squared;
    col   = idx;

    %% Row 1: 3D plot equivalent to MathieuTransitionMatrix.m
    %
    % MathieuTransitionMatrix.m uses:
    %
    % plot3(tspan, Pos2, Gesch2)
    %
    % with:
    % Pos2   = y2(2,:) = phi_dot
    % Gesch2 = y2(1,:) = phi
    %
    % Thus:
    % x-axis = psi
    % y-axis = phi_dot
    % z-axis = phi

    subplot(3,2,col);

    plot3( ...
        sol_struct(idx).tspan, ...
        sol_struct(idx).phidot_mon, ...
        sol_struct(idx).phi_mon, ...
        'b', ...
        'LineWidth', 1.4);

    grid on;
    box on;
    axis tight;

    % Equivalent axis limits used in MathieuTransitionMatrix.m
    ylim([-4.5 5.1]);
    zlim([-2.0 2.0]);

    daspect([1 1 1]);

    xlabel('$\psi \; [rad]$', ...
        'Interpreter', 'latex');

    ylabel('$\dot{\phi} \; [-]$', ...
        'Interpreter', 'latex');

    zlabel('$\phi \; [-]$', ...
        'Interpreter', 'latex');

    title(sprintf( ...
        'Second transition-matrix column, $\\nu_0^2=\\nu_C^2=%.2f$', ...
        nu_02), ...
        'Interpreter', 'latex');

    % Match MathieuTransitionMatrix.m
    view([40 15]);

    % Put the coordinate axes through the origin
    ax = gca;
    ax.XRuler.FirstCrossoverValue  = 0;
    ax.YRuler.FirstCrossoverValue  = 0;
    ax.ZRuler.FirstCrossoverValue  = 0;
    ax.ZRuler.SecondCrossoverValue = 0;
    ax.XRuler.SecondCrossoverValue = 0;
    ax.YRuler.SecondCrossoverValue = 0;

    % Match the psi tick labels
    ax.XTick = 0:pi/2:2*pi;
    ax.XTickLabel = {'0','\pi/2','\pi','3\pi/2','2\pi'};

    %% Row 2: Phase-plane plot

    subplot(3,2,col+2);

    plot( ...
        sol_struct(idx).phi_pp, ...
        sol_struct(idx).phidot_pp, ...
        'b', ...
        'LineWidth', 1.6);

    hold on;

    plot( ...
        sol_struct(idx).phi_pp(1), ...
        sol_struct(idx).phidot_pp(1), ...
        'go', ...
        'MarkerFaceColor', 'g', ...
        'MarkerSize', 8);

    plot( ...
        sol_struct(idx).phi_pp(end), ...
        sol_struct(idx).phidot_pp(end), ...
        'ro', ...
        'MarkerFaceColor', 'r', ...
        'MarkerSize', 8);

    xline(0, 'k--', 'LineWidth', 1.0);
    yline(0, 'k--', 'LineWidth', 1.0);

    hold off;

    axis equal;
    grid on;
    box on;

    xlabel('$\phi \; [-]$', ...
        'Interpreter', 'latex');

    ylabel('$\dot{\phi} \; [-]$', ...
        'Interpreter', 'latex');

    title(sprintf( ...
        'Phase plane, $\\nu_0^2=\\nu_C^2=%.2f$', ...
        nu_02), ...
        'Interpreter', 'latex');

    %% Row 3: Unwrapped angle

    subplot(3,2,col+4);

    plot( ...
        sol_struct(idx).psi_sol, ...
        sol_struct(idx).theta, ...
        'm', ...
        'LineWidth', 1.6);

    grid on;
    box on;

    xlabel('$\psi \; [rad]$', ...
        'Interpreter', 'latex');

    ylabel('$\theta(\psi)$', ...
        'Interpreter', 'latex');

    title(sprintf( ...
        '$\\theta=\\mathrm{unwrap}(\\mathrm{atan2}(\\dot{\\phi},\\phi))$, $\\nu_0^2=%.2f$', ...
        nu_02), ...
        'Interpreter', 'latex');

    ax = gca;
    ax.XTick = 0:pi/2:2*pi;
    ax.XTickLabel = {'0','\pi/2','\pi','3\pi/2','2\pi'};

end

sgtitle( ...
    'Combined winding-number visualization: $\nu_0^2=\nu_C^2=0.25$ and $5$', ...
    'Interpreter', 'latex', ...
    'FontSize', 14);

%% Figure saving is disabled as requested.
% If you later want to save, uncomment and use:
%
% fDir = 'figureFolder';
% if ~isfolder(fDir)
%     mkdir(fDir);
% end
% fileName = fullfile(fDir, 'CombinedFigure_nu0dot25_5.svg');
% exportgraphics(fig, fileName, 'ContentType', 'image');  % 'image' avoids the warning

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mathieu equation
%
% phi'' + 2D phi' + [nu_0^2 + nu_C^2 cos(psi)] phi = 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dxdpsi = MathieuDGL(psi, x, D, nu_02, nu_C2)

    K_psi = nu_02 + nu_C2*cos(psi);

    dxdpsi = [
        x(2);
        -2*D*x(2) - K_psi*x(1)
    ];

end