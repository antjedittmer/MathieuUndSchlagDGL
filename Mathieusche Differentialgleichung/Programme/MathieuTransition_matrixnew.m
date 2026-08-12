%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3D trajectory visualization of winding number
% for nu_c = nu_0 = [0.25, 5, 8]
%
% Uses external file:
%   MathieuDGL.m
%
% State convention from project:
%   x(1) = phi_dot
%   x(2) = phi
%
% 3D plot:
%   x-axis = psi
%   y-axis = phi
%   z-axis = phi_dot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
close all;
clear;

%% Configuration
isInterpTex = 0;   % 1: tex, 0: latex

if isInterpTex == 1 || verLessThan('matlab','9.8')
    isInterpTex = 1;
    strInterp = 'tex';
else
    isInterpTex = 0;
    strInterp = 'latex';
end

%% Time settings
t0 = 0.0;
T  = 2*pi;
tspan = t0:0.001:T;

%% Initial conditions (basis vectors)
Diagonal = eye(2);

%% Parameters for task 1
nu_vals = [0.25, 5, 8];
D = 0.1;   % keep fixed unless your supervisor specifies another D

%% Labels
strP.phi    = '$\phi \; [-]$';
strP.phiDot = '$\dot{\phi} \; [-]$';
strP.psi    = '$\psi \; [rad]$';

if isInterpTex == 1
    strP.phi    = '\phi [-]';
    strP.phiDot = '\phi_dot [-]';
    strP.psi    = '\psi [rad]';
end

%% Figure folder
fDir = 'figureFolder';
if ~isfolder(fDir)
    mkdir(fDir);
end

fDirTask1 = fullfile(fDir,'Winding3D');
if ~isfolder(fDirTask1)
    mkdir(fDirTask1);
end

%% Loop over nu values
for idx = 1:length(nu_vals)

    nu_02 = nu_vals(idx);
    nu_C2 = nu_02;

    %% Solve ODE for first basis vector e1 = [1;0]
    sol1 = ode45(@(psi,x) MathieuDGL(psi,x,D,nu_02,nu_C2), [t0,T], Diagonal(:,1));
    y1   = deval(sol1,tspan);

    phiDot1 = y1(1,:);
    phi1    = y1(2,:);

    %% Solve ODE for second basis vector e2 = [0;1]
    sol2 = ode45(@(psi,x) MathieuDGL(psi,x,D,nu_02,nu_C2), [t0,T], Diagonal(:,2));
    y2   = deval(sol2,tspan);

    phiDot2 = y2(1,:);
    phi2    = y2(2,:);

    %% ============================
    % Figure 1: e1 = [1;0]
    % ============================
    fig1 = figure('Color','w','Name',sprintf('nu = %.2f, e1',nu_02), ...
                  'Position',[100 100 900 700]);
    hold on;
    grid on;
    box on;

    plot3(tspan, phi1, phiDot1, 'b', 'LineWidth', 2.0);

    % axis of winding
    plot3([0 T], [0 0], [0 0], 'r--', 'LineWidth', 1.5);

    % start/end markers
    plot3(tspan(1),   phi1(1),   phiDot1(1),   'go', 'MarkerFaceColor','g', 'MarkerSize',8);
    plot3(tspan(end), phi1(end), phiDot1(end), 'ro', 'MarkerFaceColor','r', 'MarkerSize',8);

    xlabel(strP.psi, 'Interpreter', strInterp);
    ylabel(strP.phi, 'Interpreter', strInterp);
    zlabel(strP.phiDot, 'Interpreter', strInterp);

    title(sprintf('3D trajectory for $e_1=[1;0]$, $D=%.2f$, $\\nu_0=\\nu_c=%.2f$', D, nu_02), ...
          'Interpreter', strInterp);

    ax = gca;
    ax.XTick = 0:pi/2:2*pi;
    ax.XTickLabel = {'0','\pi/2','\pi','3\pi/2','2\pi'};
    ax.FontSize = 12;

    axis tight;
    daspect([1.8 1 1]);

    if abs(nu_02 - 0.25) < 1e-12
        view(40,25);
    elseif abs(nu_02 - 5) < 1e-12
        view(18,18);
        text(pi, 0.8*max(phi1), 0.8*max(phiDot1), 'expected winding: m \approx 2', ...
             'FontSize', 11, 'Interpreter', 'tex');
    else
        view(15,18);
    end

    legend({'trajectory','winding axis','start','end'}, 'Location','best');

    baseName1 = sprintf('Task1_3D_e1_D%0.2f_nu%0.2f', D, nu_02);
    baseName1 = strrep(baseName1, '.', 'dot');
    fileName1 = fullfile(fDirTask1, [baseName1 '.png']);
    exportgraphics(fig1, fileName1, 'Resolution', 300);

    %% ============================
    % Figure 2: e2 = [0;1]
    % ============================
    fig2 = figure('Color','w','Name',sprintf('nu = %.2f, e2',nu_02), ...
                  'Position',[150 120 900 700]);
    hold on;
    grid on;
    box on;

    plot3(tspan, phi2, phiDot2, 'm', 'LineWidth', 2.0);

    % axis of winding
    plot3([0 T], [0 0], [0 0], 'r--', 'LineWidth', 1.5);

    % start/end markers
    plot3(tspan(1),   phi2(1),   phiDot2(1),   'go', 'MarkerFaceColor','g', 'MarkerSize',8);
    plot3(tspan(end), phi2(end), phiDot2(end), 'ro', 'MarkerFaceColor','r', 'MarkerSize',8);

    xlabel(strP.psi, 'Interpreter', strInterp);
    ylabel(strP.phi, 'Interpreter', strInterp);
    zlabel(strP.phiDot, 'Interpreter', strInterp);

    title(sprintf('3D trajectory for $e_2=[0;1]$, $D=%.2f$, $\\nu_0=\\nu_c=%.2f$', D, nu_02), ...
          'Interpreter', strInterp);

    ax = gca;
    ax.XTick = 0:pi/2:2*pi;
    ax.XTickLabel = {'0','\pi/2','\pi','3\pi/2','2\pi'};
    ax.FontSize = 12;

    axis tight;
    daspect([1.8 1 1]);

    if abs(nu_02 - 0.25) < 1e-12
        view(40,25);
    elseif abs(nu_02 - 5) < 1e-12
        view(18,18);
    else
        view(15,18);
    end

    legend({'trajectory','winding axis','start','end'}, 'Location','best');

    baseName2 = sprintf('Task1_3D_e2_D%0.2f_nu%0.2f', D, nu_02);
    baseName2 = strrep(baseName2, '.', 'dot');
    fileName2 = fullfile(fDirTask1, [baseName2 '.png']);
    exportgraphics(fig2, fileName2, 'Resolution', 300);

    %% ============================
    % Figure 3: comparison e1 and e2
    % ============================
    fig3 = figure('Color','w','Name',sprintf('nu = %.2f, comparison',nu_02), ...
                  'Position',[200 140 950 720]);
    hold on;
    grid on;
    box on;

    plot3(tspan, phi1, phiDot1, 'b', 'LineWidth', 2.0);
    plot3(tspan, phi2, phiDot2, 'm', 'LineWidth', 2.0);
    plot3([0 T], [0 0], [0 0], 'r--', 'LineWidth', 1.5);

    xlabel(strP.psi, 'Interpreter', strInterp);
    ylabel(strP.phi, 'Interpreter', strInterp);
    zlabel(strP.phiDot, 'Interpreter', strInterp);

    title(sprintf('3D comparison, $D=%.2f$, $\\nu_0=\\nu_c=%.2f$', D, nu_02), ...
          'Interpreter', strInterp);

    ax = gca;
    ax.XTick = 0:pi/2:2*pi;
    ax.XTickLabel = {'0','\pi/2','\pi','3\pi/2','2\pi'};
    ax.FontSize = 12;

    axis tight;
    daspect([1.8 1 1]);

    if abs(nu_02 - 0.25) < 1e-12
        view(40,25);
    elseif abs(nu_02 - 5) < 1e-12
        view(18,18);
    else
        view(15,18);
    end

    legend({'e_1=[1;0]','e_2=[0;1]','winding axis'}, 'Location','best');

    baseName3 = sprintf('Task1_3D_compare_D%0.2f_nu%0.2f', D, nu_02);
    baseName3 = strrep(baseName3, '.', 'dot');
    fileName3 = fullfile(fDirTask1, [baseName3 '.png']);
    exportgraphics(fig3, fileName3, 'Resolution', 300);

end

disp('3D winding plots created.');