clc;
clearvars;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Four-bladed rotor: parameters from the slide
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Auswahl = 3;
Blatt   = 4;
AnzGl   = 2*Blatt;

gamma = 8.5;

d4 = 0.116;
d3 = 0.187;
d2 = 0.338;

ebeta = 0.15;
nu0   = 1.117;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Advance ratio
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SW    = 0.01;
MuMin = 0;
MuMax = 1;
loadData = 0
mu_paramVec = MuMin:SW:MuMax;
nMu = length(mu_paramVec);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ORIGINAL integration interval and solver settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t0 = 0;
T  = 2*pi;

Diagonal = diag(ones(AnzGl,1));

options2 = odeset( ...
    'RelTol',1e-10, ...
    'AbsTol',1e-12, ...
    'MaxStep',1e-3);

assert(exist('SchlagDGL','file') ~= 0, ...
    'Keep your original SchlagDGL.m in the same folder.');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Storage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CharMult = complex(zeros(nMu,AnzGl));

CharExRe    = zeros(nMu,AnzGl);
CharExImRaw = zeros(nMu,AnzGl);

sigmaExactUpper = zeros(nMu,1);
sigmaExactLower = zeros(nMu,1);

nuExact = zeros(nMu,1);

cond_A = zeros(nMu,1);

% Constant approximation
lambdaApprox = complex(zeros(nMu,6));
lambdaDiff   = complex(zeros(nMu,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialisation of approximation-mode tracking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% State order in your SchlagDGL.m:
%
% [beta0'; betac'; betas'; betad';
%  beta0;  betac;  betas;  betad]

idxCC   = [1 2 3 5 6 7];
idxDiff = [4 8];

[~,A0] = SchlagDGL( ...
    0,zeros(AnzGl,1), ...
    gamma,d2,d3,d4,0,ebeta,nu0,Blatt);

[V0,L0] = eig(A0(idxCC,idxCC));

lambda0 = diag(L0);

[~,order] = sort(imag(lambda0));

previousLambda = lambda0(order);
previousModes  = normalizeColumns(V0(:,order));

% Initial ordering:
% 1: negative high-frequency mode
% 2: negative collective mode
% 3: negative low-frequency mode
% 4: positive low-frequency mode
% 5: positive collective mode
% 6: positive high-frequency mode

idxHigh = [1 6];
idxColl = [2 5];
idxLow  = [3 4];

% One-to-one eigenvalue matching without additional toolboxes
P = perms(1:6);

rowIndices = repmat(1:6,size(P,1),1);

assignmentIndices = sub2ind( ...
    [6 6],rowIndices,P);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for idx = 1:nMu

    mu_param = mu_paramVec(idx);

    Monodromie = zeros(AnzGl);

    %==============================================================
    % ORIGINAL ODE integration: unchanged
    %==============================================================

    for k = 1:AnzGl

        sol = ode45( ...
            @(psi,x) SchlagDGL( ...
                psi,x,gamma,d2,d3,d4, ...
                mu_param,ebeta,nu0,Blatt), ...
            [t0,T],Diagonal(:,k),options2);

        Monodromie(:,k) = deval(sol,T);

    end

    cond_A(idx) = cond(Monodromie);

    %==============================================================
    % Multipliers and exponent post-processing
    %==============================================================

    charMult = eig(Monodromie);

    if any(~isfinite(charMult)) || ...
            any(abs(charMult) < 1e-14)

        error('Unresolved multiplier at mu = %.4f.', ...
            mu_param);

    end

    CharMult(idx,:) = charMult.';

    % Actual real parts from the multipliers
    sigma = log(abs(charMult))/T;

    CharExRe(idx,:) = sigma.';

    sigmaExactUpper(idx) = max(sigma);
    sigmaExactLower(idx) = min(sigma);

    % Principal imaginary parts:
    % use angle, not atan(imag(rho)/real(rho)).
    nuPrincipal = angle(charMult)/T;

    CharExImRaw(idx,:) = nuPrincipal.';

    %--------------------------------------------------------------
    % Slide frequency representation
    %
    % For THIS model, parameter set and range 0 <= mu <= 1,
    % full-revolution multipliers give the common frequency
    % deviation from 1/rev.
    %
    % No n, no nAdd, no automatic half-integer increments.
    %--------------------------------------------------------------

    nuExact(idx) = 1 + max(abs(nuPrincipal));

    %==============================================================
    % Constant approximation: separate additional calculation
    %==============================================================

    Naverage = 64;
    Amean = zeros(AnzGl);

    for j = 0:Naverage-1

        psi = j*T/Naverage;

        [~,Atemp] = SchlagDGL( ...
            psi,zeros(AnzGl,1), ...
            gamma,d2,d3,d4, ...
            mu_param,ebeta,nu0,Blatt);

        Amean = Amean + Atemp/Naverage;

    end

    %--------------------------------------------------------------
    % Collective/cyclic modes
    %--------------------------------------------------------------

    [V,L] = eig(Amean(idxCC,idxCC));

    currentLambda = diag(L);
    V = normalizeColumns(V);

    overlap = min(1,abs(previousModes'*V).^2);

    cost = ...
        abs(previousLambda-currentLambda.').^2 ...
        + 0.02*(1-overlap);

    totalCost = sum(cost(assignmentIndices),2);

    [~,best] = min(totalCost);

    chosen = P(best,:).';

    currentLambda = currentLambda(chosen);
    V = V(:,chosen);

    lambdaApprox(idx,:) = currentLambda.';

    previousLambda = currentLambda;
    previousModes  = V;

    %--------------------------------------------------------------
    % Differential mode
    %--------------------------------------------------------------

    ld = eig(Amean(idxDiff,idxDiff));

    [~,orderDiff] = sort(imag(ld));

    lambdaDiff(idx,:) = ld(orderDiff).';

    %--------------------------------------------------------------
    % Progress
    %--------------------------------------------------------------

    if mod(idx-1,10) == 0 || idx == nMu

        fprintf('mu = %.2f (%d/%d)\n', ...
            mu_param,idx,nMu);

    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Approximation curves
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sigmaColl = real(lambdaApprox(:,idxColl));
sigmaHigh = real(lambdaApprox(:,idxHigh));
sigmaLow  = real(lambdaApprox(:,idxLow));

sigmaDiff = real(lambdaDiff(:,2));

% Display shifts explicitly stated in the slide
nuColl = abs(imag(lambdaApprox(:,idxColl)));

nuHigh = abs(imag(lambdaApprox(:,idxHigh))) - 1;

nuLow = abs(imag(lambdaApprox(:,idxLow))) + 1;

nuDiff = abs(imag(lambdaDiff(:,2)));

muPlot = mu_paramVec(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure: only the two plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cExact = [0.00 0.00 0.55];
cColl  = [0.00 0.40 0.00];
cHigh  = [0.90 0.20 0.05];
cLow   = [0.70 0.00 0.60];
cDiff  = [0.95 0.60 0.00];

fig = figure( ...
    'Color','w', ...
    'Name','Real and imaginary characteristic exponents', ...
    'NumberTitle','off', ...
    'Position',[150 100 850 720]);

tl = tiledlayout( ...
    2,1, ...
    'TileSpacing','compact', ...
    'Padding','compact');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Real parts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ax1 = nexttile(tl);

hold(ax1,'on');
box(ax1,'on');

plot(ax1,muPlot,sigmaColl, ...
    '--','Color',cColl,'LineWidth',1.2);

plot(ax1,muPlot,sigmaHigh, ...
    '-.','Color',cHigh,'LineWidth',1.2);

plot(ax1,muPlot,sigmaLow, ...
    '--','Color',cLow,'LineWidth',1.2);

plot(ax1,muPlot,sigmaDiff, ...
    ':','Color',cDiff,'LineWidth',1.5);

plot(ax1,muPlot,sigmaExactUpper, ...
    '-','Color',cExact,'LineWidth',1.8);

plot(ax1,muPlot,sigmaExactLower, ...
    '-','Color',cExact,'LineWidth',1.8);

ylabel(ax1,'Real part \sigma [-]', ...
    'Interpreter','tex', ...
    'FontSize',15, ...
    'FontWeight','bold');

set(ax1, ...
    'XLim',[0 1], ...
    'YLim',[-0.6 -0.1], ...
    'XTick',0:0.2:1, ...
    'YTick',-0.6:0.1:-0.1, ...
    'XTickLabel',[], ...
    'FontName','Arial', ...
    'FontSize',12, ...
    'LineWidth',1.1, ...
    'TickDir','out');

grid(ax1,'off');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Imaginary parts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ax2 = nexttile(tl);

hold(ax2,'on');
box(ax2,'on');

plot(ax2,muPlot,nuColl, ...
    '--','Color',cColl,'LineWidth',1.2);

plot(ax2,muPlot,nuHigh, ...
    '-.','Color',cHigh,'LineWidth',1.2);

plot(ax2,muPlot,nuLow, ...
    '--','Color',cLow,'LineWidth',1.2);

plot(ax2,muPlot,nuDiff, ...
    ':','Color',cDiff,'LineWidth',1.5);

plot(ax2,muPlot,nuExact, ...
    '-','Color',cExact,'LineWidth',1.8);

xlabel(ax2,'Advance ratio \mu [-]', ...
    'Interpreter','tex', ...
    'FontSize',15, ...
    'FontWeight','bold');

ylabel(ax2,'Imaginary part \nu [-]', ...
    'Interpreter','tex', ...
    'FontSize',15, ...
    'FontWeight','bold');

set(ax2, ...
    'XLim',[0 1], ...
    'YLim',[0.98 1.10], ...
    'XTick',0:0.2:1, ...
    'YTick',0.98:0.02:1.10, ...
    'FontName','Arial', ...
    'FontSize',12, ...
    'LineWidth',1.1, ...
    'TickDir','out');

grid(ax2,'off');

linkaxes([ax1 ax2],'x');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Legend
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hExact = plot(ax2,NaN,NaN, ...
    '-','Color',cExact,'LineWidth',1.8);

hColl = plot(ax2,NaN,NaN, ...
    '--','Color',cColl,'LineWidth',1.2);

hHigh = plot(ax2,NaN,NaN, ...
    '-.','Color',cHigh,'LineWidth',1.2);

hLow = plot(ax2,NaN,NaN, ...
    '--','Color',cLow,'LineWidth',1.2);

hDiff = plot(ax2,NaN,NaN, ...
    ':','Color',cDiff,'LineWidth',1.5);

lgd = legend(ax2, ...
    [hExact hColl hHigh hLow hDiff], ...
    {'Floquet theory', ...
     'Collective eigenvalue', ...
     'High frequency', ...
     'Low frequency', ...
     'Differential eigenvalue'}, ...
    'Orientation','horizontal', ...
    'FontSize',10);

lgd.Layout.Tile = 'south';

% Saving Figure
exportgraphics(fig, ...
        'RealImaginaryParts_FourBladeRotor.png', ...
        'Resolution',300);

    save('Workspace_FourBladeRotor.mat', ...
        'gamma','d2','d3','d4','ebeta','nu0', ...
        'mu_paramVec','CharMult', ...
        'CharExRe','CharExImRaw', ...
        'sigmaExactUpper','sigmaExactLower','nuExact', ...
        'lambdaApprox','lambdaDiff');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Helper used only for approximation-mode matching
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function V = normalizeColumns(V)

    columnNorm = sqrt(sum(abs(V).^2,1));

    V = bsxfun(@rdivide,V,max(columnNorm,eps));

end