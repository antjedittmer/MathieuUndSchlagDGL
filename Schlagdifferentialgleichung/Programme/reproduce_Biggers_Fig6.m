% Floquet reproduction of Biggers (1974), Figure 6.
%
% Panel (a): nu = 1.1
% Panel (b): nu = 1.0
%
% White: complex-conjugate Floquet multipliers, Delta < 0.
% Grey: real Floquet multipliers, Delta >= 0.
% Black curve: numerical Floquet boundary, Delta = 0.
%
% The grey region is a real-multiplier/critical region. It is not, by
% itself, an instability classification.
%
% State vector:
%   x(1) = beta
%   x(2) = beta_dot
%
% Single-blade equation:
%
%   beta'' + c(psi)*beta' + k(psi)*beta = 0
%
% This is Biggers' single-blade flapping equation (1), written as
%
%   beta'' + gamma/8*(1 + (4*mu/3)*sin(psi))*beta'
%       + [nu^2 + gamma/8*((4*mu/3)*cos(psi)
%       + mu^2*sin(2*psi))]*beta = 0.
%
% Parameter correspondence with Biggers:
%   mu  = p       advance ratio
%   gamma = gamma blade Lock number
%   nu = v        flapping natural frequency
%   psi = Psi     rotor azimuth
%
% In Biggers' notation, primes denote differentiation with respect to
% azimuth Psi. Equation reference: Biggers, Eq. (1).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clearvars;
close all;

%% Parameters

nuList = [1.1, 1.0];

% Horizontal axis: advance ratio mu = p in Biggers.
% Vertical axis: Lock number gamma.
muVec    = 0:0.01:0.50;
gammaVec = 0:0.25:24.0;

% Faster preliminary grid, if required:
% muVec    = 0:0.02:0.50;
% gammaVec = 0:0.50:24.0;

nMu    = numel(muVec);
nGamma = numel(gammaVec);

[MU,GAMMA] = meshgrid(muVec,gammaVec);

% One rotor revolution in azimuth:
% Biggers integrates over Psi = 0,...,2*pi in the Floquet calculation.
t0 = 0;
T  = 2*pi;

% ODE solver settings.
options = odeset( ...
    'RelTol',1e-8, ...
    'AbsTol',1e-10, ...
    'MaxStep',0.01);

%% Figure setup

fig = figure( ...
    'Color','w', ...
    'Name','Biggers Figure 6 reproduction', ...
    'Position',[100 100 1200 600]);

tl = tiledlayout(fig,1,2, ...
    'TileSpacing','compact', ...
    'Padding','compact');

%% Floquet calculation

for iNu = 1:numel(nuList)

    nu = nuList(iNu);

    discMap = zeros(nGamma,nMu);

    fprintf('\n====================================================\n');
    fprintf('Biggers Fig. 6(%c), nu = %.1f\n', ...
        char('a'+iNu-1),nu);
    fprintf('====================================================\n');

    tic;

    for iGamma = 1:nGamma

        gamma = gammaVec(iGamma);

        for iMu = 1:nMu

            mu = muVec(iMu);

            % Monodromy matrix from two independent initial conditions.
            %
            % This follows Biggers' Floquet-theory procedure: integrate
            % the single-blade equation (1) over one azimuthal period
            % for each independent initial condition to obtain the state
            % transition/monodromy matrix.

            Monodromy = zeros(2,2);
            I2 = eye(2);

            for k = 1:2

                sol = ode45( ...
                    @(psi,x) BiggersFlapwiseODE( ...
                        psi,x,gamma,mu,nu), ...
                    [t0,T], ...
                    I2(:,k), ...
                    options);

                Monodromy(:,k) = deval(sol,T);

            end

            % Floquet discriminant.
            %
            % For the 2-by-2 monodromy matrix M, the characteristic
            % equation of the Floquet multipliers is
            %
            %   lambda^2 - trace(M)*lambda + det(M) = 0.
            %
            % Therefore the multiplier discriminant is
            %
            %   Delta = trace(M)^2 - 4*det(M).
            %
            % Delta < 0 gives a complex-conjugate multiplier pair.
            % Delta >= 0 gives real or repeated multipliers.
            %
            % Equation reference: Floquet transition-matrix discussion
            % following Biggers' Eq. (1); Figure 6 is constructed from
            % these Floquet results.

            trM  = trace(Monodromy);
            detM = det(Monodromy);

            Delta = trM^2 - 4*detM;

            discMap(iGamma,iMu) = real(Delta);

            % Liouville determinant check.
            %
            % The coefficient of beta_dot in Biggers' Eq. (1) is
            %
            %   c(psi) = gamma/8*(1 + (4*mu/3)*sin(psi)).
            %
            % Hence
            %
            %   integral_0^(2*pi) c(psi)dpsi = gamma*pi/4,
            %
            % and Liouville's formula gives
            %
            %   det(M) = exp(-gamma*pi/4).
            %
            % This is a consistency check on the numerical monodromy
            % matrix, not an additional equation from Biggers.

            detExpected = exp(-gamma*pi/4);

            if abs(detM-detExpected) > ...
                    1e-7*max(1,abs(detExpected))

                warning( ...
                    ['Monodromy determinant check failed at ', ...
                     'gamma = %.4g, mu = %.4g.'], ...
                    gamma,mu);

            end

        end

        if mod(iGamma,10) == 0 || iGamma == nGamma
            fprintf('Completed gamma row %d of %d\n', ...
                iGamma,nGamma);
        end

    end

    fprintf('Elapsed time: %.1f seconds\n',toc);

    %% Critical-region classification

    discTolerance = 1e-8*max(1,max(abs(discMap(:))));

    % Delta < 0:
    %   complex-conjugate Floquet multipliers.
    %
    % Delta >= 0:
    %   real or repeated Floquet multipliers.
    %
    % Biggers calls the corresponding regions critical regions. They are
    % not, by themselves, an instability classification.

    criticalMap = discMap >= -discTolerance;

    fprintf('Real-multiplier points: %d of %d\n', ...
        nnz(criticalMap),numel(criticalMap));

    fprintf('Discriminant range: [%g, %g]\n', ...
        min(discMap(:)),max(discMap(:)));

    %% Plot panel

    ax = nexttile(tl,iNu);

    imagesc(ax,muVec,gammaVec,double(criticalMap));

    set(ax,'YDir','normal');

    hold(ax,'on');

    % White: Delta < 0, complex-conjugate multipliers.
    % Grey: Delta >= 0, real or repeated multipliers.

    colormap(ax,[ ...
        1.00 1.00 1.00;
        0.68 0.68 0.68]);

    caxis(ax,[0 1]);

    % Numerical Floquet boundary Delta = 0.
    %
    % This is the numerical analogue of the boundary plotted in
    % Biggers' Figure 6.

    if min(discMap(:)) <= 0 && max(discMap(:)) >= 0

        contour(ax,MU,GAMMA,discMap,[0 0], ...
            'k','LineWidth',1.2);

    end

    %% Reference cases from Biggers

    if abs(nu-1.1) < 1e-12

        % Biggers Figure 7: case a.
        yline(ax,6,'k-.','LineWidth',1.0);

        text(ax,0.02,5.15, ...
            'case a: $\nu=1.1,\ \gamma=6$', ...
            'Interpreter','latex', ...
            'FontSize',10, ...
            'BackgroundColor','w', ...
            'Margin',2);

    else

        % Biggers Figure 8: case b.
        yline(ax,6,'k-.','LineWidth',1.0);

        % Biggers Figure 9: case c.
        yline(ax,12,'k-.','LineWidth',1.0);

        text(ax,0.02,5.15, ...
            'case b: $\nu=1.0,\ \gamma=6$', ...
            'Interpreter','latex', ...
            'FontSize',10, ...
            'BackgroundColor','w', ...
            'Margin',2);

        text(ax,0.02,11.15, ...
            'case c: $\nu=1.0,\ \gamma=12$', ...
            'Interpreter','latex', ...
            'FontSize',10, ...
            'BackgroundColor','w', ...
            'Margin',2);

    end

    %% Hover repeated-root point

    % At mu = p = 0, Biggers' Eq. (1) reduces to his hover equation (4):
    %
    %   beta'' + gamma/8*beta' + nu^2*beta = 0.
    %
    % Its characteristic equation is
    %
    %   lambda^2 + gamma/8*lambda + nu^2 = 0.
    %
    % The roots coalesce when
    %
    %   (gamma/8)^2 - 4*nu^2 = 0,
    %
    % giving gamma = 16*nu.

    gammaHover = 16*nu;

    if gammaHover <= max(gammaVec)

        plot(ax,0,gammaHover, ...
            'ko', ...
            'MarkerFaceColor','k', ...
            'MarkerSize',5);

        hoverLabel = ['$\mu=0$: $\gamma=16\nu=' ...
            num2str(gammaHover,'%.1f') '$'];

        text(ax,0.012,gammaHover+0.55, ...
            hoverLabel, ...
            'Interpreter','latex', ...
            'FontSize',9);

    end

    %% Formatting

    % Construct this title without sprintf, so LaTeX commands such as \nu
    % are not interpreted as sprintf escape sequences.

    panelTitle = ['(' char('a'+iNu-1) ') $\nu=' ...
        num2str(nu,'%.1f') '$'];

    title(ax,panelTitle, ...
        'Interpreter','latex', ...
        'FontSize',15);

    xlabel(ax,'$\mu$','Interpreter','latex','FontSize',14);
    ylabel(ax,'$\gamma$','Interpreter','latex','FontSize',14);

    xlim(ax,[0 0.5]);
    ylim(ax,[0 24]);

    xticks(ax,0:0.1:0.5);
    yticks(ax,0:4:24);

    grid(ax,'on');
    box(ax,'on');

    set(ax, ...
        'FontSize',11, ...
        'TickLabelInterpreter','latex');

end

%% Figure-level title and annotation

sgtitle(tl,{ ...
    'Reproduction of Biggers (1974), Figure 6', ...
    '$\gamma$--$\mu$ plane for a single blade in rotating coordinates, based on Floquet theory'}, ...
    'Interpreter','latex', ...
    'FontSize',14);

annotation(fig, ...
    'textbox',[0.27 0.015 0.46 0.045], ...
    'String', ...
    'White: complex-conjugate multipliers ($\Delta<0$); grey: real multipliers ($\Delta\geq0$); black: $\Delta=0$', ...
    'Interpreter','latex', ...
    'EdgeColor','none', ...
    'HorizontalAlignment','center', ...
    'FontSize',10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dx = BiggersFlapwiseODE(psi,x,gamma,mu,nu)

    % State vector:
    %   x(1) = beta
    %   x(2) = beta_dot

    beta    = x(1);
    betadot = x(2);

    % Biggers' single-blade flapping equation (1):
    %
    %   beta'' + gamma/8*(1 + (4*mu/3)*sin(psi))*beta'
    %       + [nu^2 + gamma/8*((4*mu/3)*cos(psi)
    %       + mu^2*sin(2*psi))]*beta = 0.
    %
    % The following coefficient is the beta_dot coefficient:
    %
    %   c(psi) = gamma/8*(1 + (4*mu/3)*sin(psi)).
    %
    % Equation reference: Biggers, Eq. (1).

    cPsi = gamma/8 * ...
        (1 + (4*mu/3)*sin(psi));

    % The following coefficient is the beta coefficient:
    %
    %   k(psi) = nu^2 + gamma/8*((4*mu/3)*cos(psi)
    %                              + mu^2*sin(2*psi)).
    %
    % Equation reference: Biggers, Eq. (1).

    kPsi = nu^2 + gamma/8 * ...
        ((4*mu/3)*cos(psi) + mu^2*sin(2*psi));

    % First-order state-space form of Biggers' Eq. (1):
    %
    %   x_1' = x_2,
    %   x_2' = -c(psi)*x_2 - k(psi)*x_1.

    dx = [ ...
        betadot;
        -cPsi*betadot - kPsi*beta ...
    ];

end