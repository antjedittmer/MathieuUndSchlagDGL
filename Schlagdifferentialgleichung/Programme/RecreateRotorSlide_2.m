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

mu_paramVec = MuMin:SW:MuMax;
muPlot = mu_paramVec(:);

nMu = length(mu_paramVec);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Original integration interval and solver settings
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

frequencyDeviation = zeros(nMu,1);

cond_A = zeros(nMu,1);

% Averaged-system eigenvalues
lambdaApprox = complex(zeros(nMu,6));
lambdaDiff   = complex(zeros(nMu,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialisation of averaged-system mode tracking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% State ordering in your original SchlagDGL.m:
%
% x = [beta0'; betac'; betas'; betad';
%      beta0;  betac;  betas;  betad]
%
% Collective/cyclic states:
idxCC = [1 2 3 5 6 7];

% Differential states:
idxDiff = [4 8];

[~,A0] = SchlagDGL( ...
    0,zeros(AnzGl,1), ...
    gamma,d2,d3,d4,0,ebeta,nu0,Blatt);

[V0,L0] = eig(A0(idxCC,idxCC));

lambda0 = diag(L0);

[~,initialOrder] = sort(imag(lambda0));

previousLambda = lambda0(initialOrder);

previousModes = normalizeColumns( ...
    V0(:,initialOrder));

% Initial order:
%
% 1: negative high-frequency mode
% 2: negative collective mode
% 3: negative low-frequency mode
% 4: positive low-frequency mode
% 5: positive collective mode
% 6: positive high-frequency mode

idxHigh = [1 6];
idxColl = [2 5];
idxLow  = [3 4];

% All one-to-one assignments of six eigenvalues.
% No additional MATLAB toolbox is required.
permutationList = perms(1:6);

rowList = repmat( ...
    1:6,size(permutationList,1),1);

assignmentEntries = sub2ind( ...
    [6 6],rowList,permutationList);

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
    % Floquet multipliers
    %==============================================================

    charMult = eig(Monodromie);

    if any(~isfinite(charMult)) || ...
            any(abs(charMult) < 1e-14)

        error('Unresolved multiplier at mu = %.4f.', ...
            mu_param);

    end

    CharMult(idx,:) = charMult.';

    %==============================================================
    % Real parts directly from multiplier magnitudes
    %==============================================================

    sigma = log(abs(charMult))/T;

    CharExRe(idx,:) = sigma.';

    sigmaExactUpper(idx) = max(sigma);
    sigmaExactLower(idx) = min(sigma);

    %==============================================================
    % Principal imaginary parts
    %
    % No automatic n, nAdd or half-integer increments.
    %==============================================================

    nuPrincipal = angle(charMult)/T;

    CharExImRaw(idx,:) = nuPrincipal.';

    % For this identical-blade example over 0 <= mu <= 1,
    % the full-revolution phases give the common frequency
    % deviation used to reconstruct the displayed representatives.
    %
    % This is specific to the present example, not a general
    % mode-unwrapping algorithm.
    frequencyDeviation(idx) = max(abs(nuPrincipal));

    %==============================================================
    % Constant approximation: average the original matrix
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
    % Collective/cyclic averaged subsystem
    %--------------------------------------------------------------

    [V,L] = eig(Amean(idxCC,idxCC));

    currentLambda = diag(L);
    V = normalizeColumns(V);

    overlap = min(1,abs(previousModes'*V).^2);

    cost = ...
        abs(previousLambda-currentLambda.').^2 ...
        + 0.02*(1-overlap);

    totalCost = sum( ...
        cost(assignmentEntries),2);

    [~,bestPermutation] = min(totalCost);

    assignment = ...
        permutationList(bestPermutation,:).';

    currentLambda = currentLambda(assignment);
    V = V(:,assignment);

    lambdaApprox(idx,:) = currentLambda.';

    previousLambda = currentLambda;
    previousModes  = V;

    %--------------------------------------------------------------
    % Differential averaged subsystem
    %--------------------------------------------------------------

    differentialEigenvalues = eig( ...
        Amean(idxDiff,idxDiff));

    [~,orderDiff] = sort( ...
        imag(differentialEigenvalues));

    lambdaDiff(idx,:) = ...
        differentialEigenvalues(orderDiff).';

    %--------------------------------------------------------------
    % Progress
    %--------------------------------------------------------------

    if mod(idx-1,10) == 0 || idx == nMu

        fprintf('mu = %.2f (%d/%d)\n', ...
            mu_param,idx,nMu);

    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cLow  = [0.10 0.10 0.65];
cColl = [0.08 0.45 0.10];
cHigh = [0.70 0.02 0.58];

lineExact  = 1.6;
lineApprox = 1.3;

ShowModeLabels = true;
ShowSplitLabel = true;

SaveFigure  = true;
SaveResults = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Complex-plane figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig = figure( ...
    'Color','w', ...
    'Name','Characteristic exponents in the complex plane', ...
    'NumberTitle','off', ...
    'Position',[150 80 780 780]);

% Frequency bands, ordered from top to bottom
bandCentre = [2 1 0 -1 -2];

bandLimits = [
     1.99   2.09
     0.99   1.09
    -0.09   0.09
    -1.09  -0.99
    -2.09  -1.99
];

bandTicks = {
    [2.00 2.08]
    [1.00 1.08]
    [-0.08 0 0.08]
    [-1.08 -1.00]
    [-2.08 -2.00]
};

bandTickLabels = {
    {'2.00','2.08'}
    {'1.00','1.08'}
    {'-0.08','0','0.08'}
    {'-1.08','-1.00'}
    {'-2.08','-2.00'}
};

left  = 0.12;
width = 0.70;

bottom = [0.79 0.63 0.40 0.24 0.08];
height = [0.13 0.13 0.20 0.13 0.13];

bandColours = [
    cHigh
    cColl
    cLow
    cColl
    cHigh
];

ax = gobjects(5,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Draw all five frequency bands
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for b = 1:5

    ax(b) = axes( ...
        'Parent',fig, ...
        'Position',[left bottom(b) width height(b)]);

    hold(ax(b),'on');

    centre = bandCentre(b);
    colour = bandColours(b,:);

    %==============================================================
    % Constant approximation: dashed
    % Use the unshifted complex eigenvalues.
    %==============================================================

    switch b

        case {1,5}

            approximation = ...
                lambdaApprox(:,idxHigh);

        case {2,4}

            approximation = [ ...
                lambdaApprox(:,idxColl), ...
                lambdaDiff];

        case 3

            approximation = ...
                lambdaApprox(:,idxLow);

    end

    for j = 1:size(approximation,2)

        plot(ax(b), ...
            real(approximation(:,j)), ...
            imag(approximation(:,j)), ...
            '--', ...
            'Color',colour, ...
            'LineWidth',lineApprox);

    end

    %==============================================================
    % Exact Floquet representatives: solid
    %
    % Positive bands: centre + frequencyDeviation
    % Negative bands: centre - frequencyDeviation
    % Zero band: both signs
    %
    % Integer-shifted representatives are displayed to match the
    % lecture figure. They are not additional physical modes.
    %==============================================================

    if centre > 0

        frequencySigns = 1;

    elseif centre < 0

        frequencySigns = -1;

    else

        frequencySigns = [-1 1];

    end

    for sgn = frequencySigns

        yExact = centre + sgn*frequencyDeviation;

        plot(ax(b), ...
            sigmaExactUpper,yExact, ...
            '-', ...
            'Color',colour, ...
            'LineWidth',lineExact);

        plot(ax(b), ...
            sigmaExactLower,yExact, ...
            '-', ...
            'Color',colour, ...
            'LineWidth',lineExact);

        %----------------------------------------------------------
        % Direction arrows from computed trajectory points
        %----------------------------------------------------------

        k1 = max(1,nMu-10);
        k2 = nMu;

        dampingBranches = [ ...
            sigmaExactUpper, ...
            sigmaExactLower];

        for branch = 1:2

            x1 = dampingBranches(k1,branch);
            x2 = dampingBranches(k2,branch);

            y1 = yExact(k1);
            y2 = yExact(k2);

            insideBand = ...
                y1 >= bandLimits(b,1) && ...
                y1 <= bandLimits(b,2) && ...
                y2 >= bandLimits(b,1) && ...
                y2 <= bandLimits(b,2);

            if insideBand && hypot(x2-x1,y2-y1) > 1e-12

                quiver(ax(b), ...
                    x1,y1,x2-x1,y2-y1,0, ...
                    'Color',colour, ...
                    'LineWidth',1.3, ...
                    'MaxHeadSize',0.8, ...
                    'HandleVisibility','off');

            end

        end

    end

    %==============================================================
    % Axis styling
    %==============================================================

    set(ax(b), ...
        'XLim',[-0.60 -0.10], ...
        'YLim',bandLimits(b,:), ...
        'XTick',-0.6:0.1:-0.1, ...
        'YTick',bandTicks{b}, ...
        'YTickLabel',bandTickLabels{b}, ...
        'YAxisLocation','right', ...
        'FontName','Arial', ...
        'FontSize',12, ...
        'LineWidth',1.1, ...
        'TickDir','out', ...
        'Box','off', ...
        'Color',[0.96 0.96 0.96]);

    grid(ax(b),'off');

    if b ~= 5

        set(ax(b),'XColor','none');

    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Zero-frequency reference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(ax(3), ...
    [-0.60 -0.10],[0 0], ...
    '-', ...
    'Color',[0.35 0.35 0.35], ...
    'LineWidth',0.8);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Axis labels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xlabel(ax(5), ...
    'Real part \sigma [-]', ...
    'Interpreter','tex', ...
    'FontName','Arial', ...
    'FontSize',16, ...
    'FontWeight','bold');

ylabel(ax(3), ...
    'Imaginary part \nu [-]', ...
    'Interpreter','tex', ...
    'FontName','Arial', ...
    'FontSize',16, ...
    'FontWeight','bold');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Break marks on the frequency axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for b = 1:4

    xBreak = left + width;

    gapTop = bottom(b);

    gapBottom = ...
        bottom(b+1) + height(b+1);

    yBreak = (gapTop+gapBottom)/2;

    annotation(fig,'line', ...
        [xBreak-0.007 xBreak+0.007], ...
        [yBreak+0.010 yBreak], ...
        'Color',[0.2 0.2 0.2], ...
        'LineWidth',1.2);

    annotation(fig,'line', ...
        [xBreak-0.007 xBreak+0.007], ...
        [yBreak yBreak-0.010], ...
        'Color',[0.2 0.2 0.2], ...
        'LineWidth',1.2);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Labels inside the graph
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ShowModeLabels

    labelColour = [0.35 0.35 0.35];

    text(ax(1),-0.585,2.045, ...
        {'High frequency','eigenvalue'}, ...
        'FontSize',10, ...
        'Color',labelColour);

    text(ax(2),-0.585,1.045, ...
        {'Collective and','difference','eigenvalue'}, ...
        'FontSize',10, ...
        'Color',labelColour);

    text(ax(3),-0.585,0.055, ...
        {'Low frequency','eigenvalue'}, ...
        'FontSize',10, ...
        'Color',labelColour);

    text(ax(4),-0.445,-1.065, ...
        {'Collective and','difference','eigenvalue'}, ...
        'FontSize',10, ...
        'Color',labelColour);

    text(ax(5),-0.445,-2.065, ...
        {'High frequency','eigenvalue'}, ...
        'FontSize',10, ...
        'Color',labelColour);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Splitting-point annotation estimated from calculated results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ShowSplitLabel

    dampingSeparation = abs( ...
        sigmaExactUpper-sigmaExactLower);

    splitIndex = find( ...
        dampingSeparation > 1e-4,1,'first');

    if ~isempty(splitIndex) && splitIndex > 1

        labelIndex = splitIndex-1;

        splitSigma = 0.5*( ...
            sigmaExactUpper(labelIndex) ...
            + sigmaExactLower(labelIndex));

        text(ax(4), ...
            splitSigma+0.008,-0.997, ...
            sprintf('\\mu \\approx %.2f',muPlot(labelIndex)), ...
            'Interpreter','tex', ...
            'FontSize',10, ...
            'Color',[0.4 0.4 0.4], ...
            'VerticalAlignment','top');

    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if SaveFigure

    exportgraphics(fig, ...
        'ComplexPlane_FourBladeRotor.png', ...
        'Resolution',300);

    savefig(fig, ...
        'ComplexPlane_FourBladeRotor.fig');

    fprintf('\nFigure saved as PNG and MATLAB FIG.\n');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Optional saving of calculated results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if SaveResults

    save('Workspace_ComplexPlane_FourBladeRotor.mat', ...
        'Auswahl','Blatt', ...
        'gamma','d2','d3','d4','ebeta','nu0', ...
        'mu_paramVec','T', ...
        'CharMult','CharExRe','CharExImRaw', ...
        'sigmaExactUpper','sigmaExactLower', ...
        'frequencyDeviation', ...
        'lambdaApprox','lambdaDiff','cond_A');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Helper used only for averaged-system eigenvector matching
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function V = normalizeColumns(V)

    columnNorm = sqrt(sum(abs(V).^2,1));

    V = bsxfun( ...
        @rdivide,V,max(columnNorm,eps));

end