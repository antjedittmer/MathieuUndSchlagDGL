clc;
clearvars;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Biggers-style critical-region map for the flapwise single-blade ODE
%
% Uses SchlagDGL.m with Blatt = 1 and parameters chosen to reproduce the
% single-blade Biggers form:
%
%   e_beta = 0, B = 1
%   d2 = B^2/4
%   d3 = B^3/3
%   d4 = B^4/4
%
% Equation represented by SchlagDGL.m for Blatt = 1:
%
% beta'' + M_beta_dot(psi)*beta' + (nu0^2 + M_beta(psi))*beta = 0.
%
% The map is calculated with Floquet theory:
% 1. Integrate both columns of the identity matrix over 0 <= psi <= 2*pi.
% 2. Form the monodromy matrix.
% 3. Compute the Floquet multipliers.
% 4. Use the discriminant to detect real/split multipliers.
%
% Figure colours:
%   white  = complex-conjugate multiplier pair (ordinary region)
%   grey   = real/split multiplier pair (critical region)
%   blue   = positive-real multiplier class: 0/rev modulo 1
%   orange = negative-real multiplier class: 1/2-rev modulo 1
%
% IMPORTANT:
% A Floquet multiplier determines Im(s) only modulo 1/rev.
% Therefore raw multiplier data cannot distinguish 0/rev from 1/rev.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Model and integration settings

Blatt = 1;

% Biggers-style single-blade parameters.
ebeta = 0;
B = 1;

d2 = B^2/4;
d3 = B^3/3;
d4 = B^4/4;

t0 = 0;
T = 2*pi;

% First use this medium grid. It is much faster than a very fine grid.
muVec = 0:0.01:0.5;
gammaVec = 0.2:0.25:24;

nMu = numel(muVec);
nGamma = numel(gammaVec);

% Use nu = 1.1 and nu = 1.0, corresponding to Biggers Fig. 6.
nuList = [1.1, 1.0];

options = odeset( ...
    'RelTol',1e-8, ...
    'AbsTol',1e-10, ...
    'MaxStep',0.01);

% Numerical tolerance for calling a multiplier real.
imagTolerance = 1e-6;

%% One common figure

figure('Color','w','Position',[100 100 1150 550]);

for iNu = 1:numel(nuList)

    nu0 = nuList(iNu);

    % regionCode:
    % 0 = ordinary: complex-conjugate multipliers
    % 1 = critical: positive-real multipliers
    % 2 = critical: negative-real multipliers
    regionCode = zeros(nGamma,nMu);

    % Discriminant of monodromy characteristic polynomial.
    % Positive: two real multipliers; negative: complex pair.
    discMap = zeros(nGamma,nMu);

    % Largest real part of characteristic exponent, diagnostic only.
    sigmaMap = zeros(nGamma,nMu);

    tic;

    for iG = 1:nGamma

        gamma = gammaVec(iG);

        for iM = 1:nMu

            mu_param = muVec(iM);

            % Build the monodromy matrix with the same independent-column
            % procedure as in BerechnungSchlagDGL.m.
            Monodromie = zeros(2,2);
            Diagonal = eye(2);

            for k = 1:2
                sol = ode45( ...
                    @(psi,x) SchlagDGL( ...
                        psi,x,gamma,d2,d3,d4, ...
                        mu_param,ebeta,nu0,Blatt), ...
                    [t0,T],Diagonal(:,k),options);

                Monodromie(:,k) = deval(sol,T);
            end

            % Floquet multipliers.
            charMult = eig(Monodromie);

            % Floquet characteristic exponents.
            charExp = log(charMult)/T;
            sigmaMap(iG,iM) = max(real(charExp));

            % Discriminant: determines whether the multiplier pair is real.
            discValue = real(trace(Monodromie)^2 - 4*det(Monodromie));
            discMap(iG,iM) = discValue;

            isRealPair = abs(imag(charMult(1))) < imagTolerance && ...
                         abs(imag(charMult(2))) < imagTolerance;

            if isRealPair
                % In a real/split region, the multipliers normally have the
                % same sign because det(Phi(T)) > 0 for this ODE.
                if mean(real(charMult)) >= 0
                    regionCode(iG,iM) = 1;
                else
                    regionCode(iG,iM) = 2;
                end
            end
        end
    end

    fprintf('nu0 = %.2f: %.1f seconds\n',nu0,toc);

    %% Plot

    subplot(1,2,iNu);
    hold on;

    % imagesc avoids warnings if all points have the same class.
    imagesc(muVec,gammaVec,regionCode);
    set(gca,'YDir','normal');

    % 0: white, 1: blue, 2: orange.
    colormap(gca,[ ...
        1.00 1.00 1.00;
        0.35 0.60 0.90;
        0.95 0.60 0.15]);

    caxis([0 2]);

    % Numerical critical-region boundaries:
    % draw only when discriminant = 0 occurs inside the grid.
    if min(discMap(:)) <= 0 && max(discMap(:)) >= 0
        contour(muVec,gammaVec,discMap,[0 0], ...
            'k','LineWidth',1.1);
    end

    % At mu = 0, the ordinary damping transition is gamma = 16*nu0.
    gammaHover = 16*nu0;
    if gammaHover >= min(gammaVec) && gammaHover <= max(gammaVec)
        plot(0,gammaHover,'ko','MarkerFaceColor','k','MarkerSize',5);
        text(0.01,gammaHover+0.55, ...
            sprintf('$\\gamma=16\\nu=%.1f$',gammaHover), ...
            'Interpreter','latex','FontSize',9);
    end

    % CASE reference lines from Biggers Figs. 7, 8, and 9.
    if abs(nu0-1.1) < 1e-12
        yline(6,'-.','Color',[0.25 0.25 0.25],'LineWidth',1);
        text(0.02,5.2,'CASE a: $\nu=1.1,\gamma=6$', ...
            'Interpreter','latex','FontSize',9);
    elseif abs(nu0-1.0) < 1e-12
        yline(6,'-.','Color',[0.25 0.25 0.25],'LineWidth',1);
        yline(12,'-.','Color',[0.25 0.25 0.25],'LineWidth',1);
        text(0.02,5.2,'CASE b: $\nu=1.0,\gamma=6$', ...
            'Interpreter','latex','FontSize',9);
        text(0.02,11.2,'CASE c: $\nu=1.0,\gamma=12$', ...
            'Interpreter','latex','FontSize',9);
    end

    xlabel('Advance ratio $\mu$','Interpreter','latex');
    ylabel('Lock number $\gamma$','Interpreter','latex');
    title(sprintf('Biggers-style critical regions: $\\nu=%.1f$',nu0), ...
        'Interpreter','latex');

    xlim([0 0.5]);
    ylim([0 24]);

    xticks(0:0.1:0.5);
    yticks(0:4:24);

    grid on;
    box on;
    set(gca,'FontSize',11,'TickLabelInterpreter','latex');

    % Diagnostic output.
    nCritical = nnz(regionCode > 0);
    fprintf('  Critical grid points: %d of %d\n', ...
        nCritical,numel(regionCode));
    fprintf('  discriminant range: [%g, %g]\n', ...
        min(discMap(:)),max(discMap(:)));
end

% A shared legend in an annotation box.
annotation('textbox',[0.30 0.01 0.42 0.06], ...
    'String', ...
    {'white: complex-conjugate multipliers (ordinary)', ...
     'blue: positive-real multipliers (0/rev modulo 1); orange: negative-real multipliers (1/2-rev modulo 1)'}, ...
    'EdgeColor','none', ...
    'HorizontalAlignment','center', ...
    'FontSize',9);

sgtitle({ ...
    'Flapwise blade ODE: Biggers-style Floquet critical-region map', ...
    'Black curves: numerical boundary where Floquet multipliers coalesce'}, ...
    'FontWeight','bold');