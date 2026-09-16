%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Combined Winding Number plots
%
% Figure 1: 3x2 (all plots)
% Figure 2: 1x2 (only 3D plots)
% Figure 3: 2x2 (only 2D plots: phase plane + angle)
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

%% Common 3D axis limits (same for both nu^2 values)
xlim3D = [0, 2*pi];
ylim3D = [-4.5, 5.1];
zlim3D = [-2.0, 2.0];

%% ------------------------------------------------------------------------
%% Figure 1: 3x2 layout (all plots, as before)
%% ------------------------------------------------------------------------
fig1 = figure( ...
    'Color', 'w', ...
    'Name', 'Combined winding-number figure (3x2)', ...
    'Position', [100 100 1300 900]);

for idx = 1:length(nu_squared_vals)
    nu_02 = sol_struct(idx).nu_squared;
    col   = idx;
    
    %% Row 1: 3D plot (second transition-matrix column)
    subplot(3,2,col);
    plot3( ...
        sol_struct(idx).tspan, ...
        sol_struct(idx).phidot_mon, ...
        sol_struct(idx).phi_mon, ...
        'b', ...
        'LineWidth', 1.4);
    grid off;
    box off;  % no box, "book style"
    
    xlim(xlim3D);
    ylim(ylim3D);
    zlim(zlim3D);
    daspect([1 1 1]);
    
    xlabel('$\psi \; [rad]$', 'Interpreter', 'latex');
    ylabel('$\dot{\phi} \; [-]$', 'Interpreter', 'latex');
    zlabel('$\phi \; [-]$', 'Interpreter', 'latex');
    title(sprintf( ...
        'Second transition-matrix column, $\\nu_0^2=\\nu_C^2=%.2f$', ...
        nu_02), ...
        'Interpreter', 'latex');
    view([40 15]);
    
    % Axes through origin
    ax = gca;
    ax.XRuler.FirstCrossoverValue  = 0;
    ax.YRuler.FirstCrossoverValue  = 0;
    ax.ZRuler.FirstCrossoverValue  = 0;
    ax.ZRuler.SecondCrossoverValue = 0;
    ax.XRuler.SecondCrossoverValue = 0;
    ax.YRuler.SecondCrossoverValue = 0;
    
    % Psi tick labels
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

    % axis tight;
    % axis padded;
    
    grid on;
    box on;
    
    xlabel('$\phi \; [-]$', 'Interpreter', 'latex');
    ylabel('$\dot{\phi} \; [-]$', 'Interpreter', 'latex');
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

    % axis tight;
    % axis padded;

    grid on;
    box on;
    
    xlabel('$\psi \; [rad]$', 'Interpreter', 'latex');
    ylabel('$\theta(\psi)$', 'Interpreter', 'latex');
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

%% ------------------------------------------------------------------------
%% Figure 2: 1x2 layout (only the two 3D plots)
%% ------------------------------------------------------------------------
fig2 = figure( ...
    'Color', 'w', ...
    'Name', '3D plots only (1x2)', ...
    'Position', [100 100 1300 600]);

for idx = 1:length(nu_squared_vals)
    nu_02 = sol_struct(idx).nu_squared;
    
    subplot(1,2,idx);
    plot3( ...
        sol_struct(idx).tspan, ...
        sol_struct(idx).phidot_mon, ...
        sol_struct(idx).phi_mon, ...
        'b', ...
        'LineWidth', 1.4);
    grid off;
    box off;
    
    xlim(xlim3D);
    ylim(ylim3D);
    zlim(zlim3D);
    daspect([1 1 1]);
    
    xlabel('$\psi \; [rad]$', 'Interpreter', 'latex');
    ylabel('$\dot{\phi} \; [-]$', 'Interpreter', 'latex');
    zlabel('$\phi \; [-]$', 'Interpreter', 'latex');
    title(sprintf( ...
        'Second transition-matrix column, $\\nu_0^2=\\nu_C^2=%.2f$', ...
        nu_02), ...
        'Interpreter', 'latex');
    view([40 15]);
    
    ax = gca;
    ax.XRuler.FirstCrossoverValue  = 0;
    ax.YRuler.FirstCrossoverValue  = 0;
    ax.ZRuler.FirstCrossoverValue  = 0;
    ax.ZRuler.SecondCrossoverValue = 0;
    ax.XRuler.SecondCrossoverValue = 0;
    ax.YRuler.SecondCrossoverValue = 0;
    ax.XTick = 0:pi/2:2*pi;
    ax.XTickLabel = {'0','\pi/2','\pi','3\pi/2','2\pi'};
end

sgtitle( ...
    '3D winding-number plots: second transition-matrix column', ...
    'Interpreter', 'latex', ...
    'FontSize', 14);

%% ------------------------------------------------------------------------
%% Figure 3: 2x2 layout (only 2D plots: phase plane + angle)
%% ------------------------------------------------------------------------
fig3 = figure( ...
    'Color', 'w', ...
    'Name', '2D plots only (2x2)', ...
    'Position', [100 100 1000 800]);

for idx = 1:length(nu_squared_vals)
    nu_02 = sol_struct(idx).nu_squared;
    
    % Phase plane
    subplot(2,2,2*idx-1);
    plot( ...
        sol_struct(idx).phi_pp, ...
        sol_struct(idx).phidot_pp, ...
        'b', 'LineWidth', 1.6);
    hold on;
    plot(sol_struct(idx).phi_pp(1), sol_struct(idx).phidot_pp(1), ...
        'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
    plot(sol_struct(idx).phi_pp(end), sol_struct(idx).phidot_pp(end), ...
        'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
    xline(0, 'k--', 'LineWidth', 1.0);
    yline(0, 'k--', 'LineWidth', 1.0);
    hold off;

    % axis tight;
    % axis padded;
    
    grid on;
    box on;
    
    xlabel('$\phi \; [-]$', 'Interpreter', 'latex');
    ylabel('$\dot{\phi} \; [-]$', 'Interpreter', 'latex');
    title(sprintf( ...
        'Phase plane, $\\nu_0^2=\\nu_C^2=%.2f$', ...
        nu_02), ...
        'Interpreter', 'latex');
    
    % Unwrapped angle
    subplot(2,2,2*idx);
    plot( ...
        sol_struct(idx).psi_sol, ...
        sol_struct(idx).theta, ...
        'm', 'LineWidth', 1.6);

    % axis tight;
    % axis padded;
    
    grid on;
    box on;
    
    xlabel('$\psi \; [rad]$', 'Interpreter', 'latex');
    ylabel('$\theta(\psi)$', 'Interpreter', 'latex');
    title(sprintf( ...
        '$\\theta=\\mathrm{unwrap}(\\mathrm{atan2}(\\dot{\\phi},\\phi))$, $\\nu_0^2=%.2f$', ...
        nu_02), ...
        'Interpreter', 'latex');
    
    ax = gca;
    ax.XTick = 0:pi/2:2*pi;
    ax.XTickLabel = {'0','\pi/2','\pi','3\pi/2','2\pi'};
end

sgtitle('2D winding-number plots: phase plane and unwrapped angle', ...
    'Interpreter', 'latex', 'FontSize', 14);

%% ------------------------------------------------------------------------
%% Save all three figures as SVG files
%% ------------------------------------------------------------------------
%% Save figures as SVG files

figureFolder = 'figureFolder';

if ~isfolder(figureFolder)
    mkdir(figureFolder);
end

saveas(fig1, fullfile(figureFolder, 'CombinedWinding_3x2.svg'));
saveas(fig2, fullfile(figureFolder, 'Winding_3D_1x2.svg'));
saveas(fig3, fullfile(figureFolder, 'Winding_2D_2x2.svg'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mathieu equation
%
% phi'' + 2D phi' + [nu_0^2 + nu_C^2 cos(psi)] phi = 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dxdpsi = MathieuDGL(psi, x, D, nu_02, nu_C2)
    % State-space form of the damped Mathieu equation
    K_psi = nu_02 + nu_C2*cos(psi);
    dxdpsi = [
        x(2);
        -2*D*x(2) - K_psi*x(1)
    ];
end