%% 
clc;
clearvars;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Rotor-Auswahl
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AuswahlInfo = {
    1, 3, '3-Blatt-Rotor, see-saw';
    2, 3, '3-Blatt-Rotor, voll gelenkig';
    3, 4, '4-Blatt-Rotor, voll gelenkig';
    4, 5, '5-Blatt-Rotor, voll gelenkig';
    5, 3, '3-Blatt-Rotor, gelenk-/lagerlos';
    6, 4, '4-Blatt-Rotor, gelenk-/lagerlos';
    7, 1, 'Einzelblattkoordinaten im rotierenden System'
    };

% Four-bladed rotor from the slide
Auswahl = 3;

Blatt = AuswahlInfo{Auswahl, 2};
AnzGl = 2 * Blatt;

fprintf('Auswahl %d: %s\n', ...
    Auswahl, AuswahlInfo{Auswahl, 3});

% The supplied matrices distinguish blade counts.
% This script does not add separate see-saw or bearingless constraints.
assert(ismember(Auswahl, [2 3 4 7]), ...
    'Use Auswahl = 2, 3, 4 or 7 with the supplied matrices.');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter aus der Vorlesungsfolie
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gamma = 8.5;

d4 = 0.116;
d3 = 0.187;
d2 = 0.338;

ebeta = 0.15;
nu0   = 1.117;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fortschrittsgrad
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SW    = 0.01;
MuMin = 0;
MuMax = 1;

mu_paramVec = (MuMin:SW:MuMax).';
nMu = numel(mu_paramVec);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Integrationsintervall
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t0 = 0;

% Periods of the matrices in the supplied SchlagDGL.m
switch Blatt

    case 1
        T = 2*pi;

    case 3
        T = 2*pi/3;

    case 4
        T = pi;

    case 5
        T = 2*pi/5;

    otherwise
        error('Unsupported blade count.');

end

branchSpacing = 2*pi/T;

options = odeset( ...
    'RelTol', 1e-10, ...
    'AbsTol', 1e-12, ...
    'MaxStep', T/100);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Systemmatrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Afun = @(psi, mu) getSystemMatrix( ...
    psi, mu, gamma, d2, d3, d4, ebeta, nu0, Blatt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialisierung bei mu = 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A0 = Afun(0, 0);

[V0, L0] = eig(A0);
lambda0 = diag(L0);

% Initial ordering by frequency
[~, initialOrder] = sort(imag(lambda0));

lambda0 = lambda0(initialOrder);
V0 = normalizeColumns(V0(:, initialOrder));

fprintf('\nInitial characteristic exponents:\n');
disp(lambda0);

previousLambda = lambda0;
previousModes  = V0;

previousMeanLambda = lambda0;
previousMeanModes  = V0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Speicher
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CharEx = complex(zeros(nMu, AnzGl));

CharExMean = complex(zeros(nMu, AnzGl));

CharMult = complex(zeros(nMu, AnzGl));

MonodromieAll = zeros(AnzGl, AnzGl, nMu);

cond_A = zeros(nMu, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Berechnung
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for idx = 1:nMu

    mu_param = mu_paramVec(idx);

    %--------------------------------------------------------------
    % Fundamental matrix
    %--------------------------------------------------------------

    Phi0 = eye(AnzGl);

    rhs = @(psi, y) fundamentalMatrixRHS( ...
        psi, y, Afun, mu_param, AnzGl);

    [~, Y] = ode45( ...
        rhs, [t0 T], Phi0(:), options);

    Monodromie = reshape( ...
        Y(end, :).', AnzGl, AnzGl);

    MonodromieAll(:, :, idx) = Monodromie;

    cond_A(idx) = cond(Monodromie);

    %--------------------------------------------------------------
    % Characteristic multipliers
    %--------------------------------------------------------------

    [V, R] = eig(Monodromie);

    charMult = diag(R);
    V = normalizeColumns(V);

    if any(~isfinite(charMult)) || ...
            any(abs(charMult) < 1e-14)

        error( ...
            'Unresolved multiplier at mu = %.4f.', ...
            mu_param);

    end

    %--------------------------------------------------------------
    % Principal characteristic exponents
    %--------------------------------------------------------------

    principal = ( ...
        log(abs(charMult)) + ...
        1i * angle(charMult)) / T;

    %--------------------------------------------------------------
    % Individual logarithm branches
    %--------------------------------------------------------------

    shifts = round( ...
        (imag(previousLambda) - imag(principal).') ...
        / branchSpacing);

    candidates = principal.' ...
        + 1i * branchSpacing * shifts;

    %--------------------------------------------------------------
    % Match modes to the previous parameter step
    %--------------------------------------------------------------

    overlap = abs(previousModes' * V).^2;
    overlap = min(overlap, 1);

    cost = abs(candidates - previousLambda).^2 ...
        + 0.02 * (1 - overlap);

    assignment = minimumAssignment(cost);

    chosenEntries = sub2ind( ...
        [AnzGl AnzGl], ...
        (1:AnzGl).', assignment);

    lambda = candidates(chosenEntries);

    CharEx(idx, :) = lambda.';

    CharMult(idx, :) = charMult(assignment).';

    previousLambda = lambda;
    previousModes  = V(:, assignment);

    %--------------------------------------------------------------
    % Constant approximation:
    % average the supplied system matrix over one period
    %--------------------------------------------------------------

    Naverage = 64;

    Amean = zeros(AnzGl);

    for k = 0:Naverage-1

        psi = k * T / Naverage;

        Amean = Amean ...
            + Afun(psi, mu_param) / Naverage;

    end

    [Vm, Lm] = eig(Amean);

    meanLambda = diag(Lm);
    Vm = normalizeColumns(Vm);

    meanOverlap = abs(previousMeanModes' * Vm).^2;
    meanOverlap = min(meanOverlap, 1);

    meanCost = ...
        abs(previousMeanLambda - meanLambda.').^2 ...
        + 0.02 * (1 - meanOverlap);

    meanAssignment = minimumAssignment(meanCost);

    meanLambda = meanLambda(meanAssignment);
    Vm = Vm(:, meanAssignment);

    CharExMean(idx, :) = meanLambda.';

    previousMeanLambda = meanLambda;
    previousMeanModes  = Vm;

    %--------------------------------------------------------------
    % Progress
    %--------------------------------------------------------------

    if mod(idx-1, 10) == 0 || idx == nMu

        fprintf('mu = %.2f   (%d/%d)\n', ...
            mu_param, idx, nMu);

    end

end

CharExRe = real(CharEx);
CharExIm = imag(CharEx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot: only the complex-plane graph from the slide
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This particular broken-axis layout is for the four-blade example.
if Blatt ~= 4

    figure('Color', 'w');

    plot( ...
        real(CharEx), imag(CharEx), ...
        'LineWidth', 1.5);

    xlabel('Real part \sigma [-]');
    ylabel('Imaginary part \nu [-]');

    grid on;
    box on;

else

    fig = figure( ...
        'Color', 'w', ...
        'Name', 'Characteristic exponent trajectories', ...
        'NumberTitle', 'off', ...
        'Position', [150 80 790 740]);

    %--------------------------------------------------------------
    % Broken frequency axis: top to bottom
    %--------------------------------------------------------------

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
        {'2.00', '2.08'}
        {'1.00', '1.08'}
        {'-0.08', '0', '0.08'}
        {'-1.08', '-1.00'}
        {'-2.08', '-2.00'}
        };

    left  = 0.13;
    width = 0.69;

    bottom = [0.79 0.63 0.40 0.24 0.08];
    height = [0.13 0.13 0.20 0.13 0.13];

    %--------------------------------------------------------------
    % Colours
    %--------------------------------------------------------------

    colours = [
        0.10  0.15  0.70
        0.10  0.50  0.12
        0.70  0.02  0.60
        ];

    initialFrequency = abs(imag(lambda0));

    group = ones(AnzGl, 1);

    group(initialFrequency > 0.5 & ...
        initialFrequency < 1.5) = 2;

    group(initialFrequency >= 1.5) = 3;

    %--------------------------------------------------------------
    % Draw curves in each frequency band
    %--------------------------------------------------------------

    ax = gobjects(5, 1);

    for b = 1:5

        ax(b) = axes( ...
            'Parent', fig, ...
            'Position', ...
            [left bottom(b) width height(b)]);

        hold(ax(b), 'on');

        % Dashed: constant approximation
        for j = 1:AnzGl

            c = colours(group(j), :);

            % Include conjugate representatives so that both
            % positive and negative frequency bands are shown.
            % These are display copies, not additional modes.
            for sgn = [-1 1]

                plot(ax(b), ...
                    real(CharExMean(:, j)), ...
                    sgn * imag(CharExMean(:, j)), ...
                    '--', ...
                    'Color', c, ...
                    'LineWidth', 1.3);

            end

        end

        % Solid: Floquet solution
        for j = 1:AnzGl

            c = colours(group(j), :);

            for sgn = [-1 1]

                plot(ax(b), ...
                    real(CharEx(:, j)), ...
                    sgn * imag(CharEx(:, j)), ...
                    '-', ...
                    'Color', c, ...
                    'LineWidth', 1.6);

            end

        end

        set(ax(b), ...
            'XLim', [-0.60 -0.10], ...
            'YLim', bandLimits(b, :), ...
            'XTick', -0.6:0.1:-0.1, ...
            'YTick', bandTicks{b}, ...
            'YTickLabel', bandTickLabels{b}, ...
            'YAxisLocation', 'right', ...
            'FontName', 'Arial', ...
            'FontSize', 13, ...
            'LineWidth', 1.0, ...
            'TickDir', 'out', ...
            'Box', 'off', ...
            'Color', [0.96 0.96 0.96]);

        grid(ax(b), 'off');

        if b ~= 5
            set(ax(b), 'XColor', 'none');
        end

    end

    %--------------------------------------------------------------
    % Zero-frequency reference
    %--------------------------------------------------------------

    plot(ax(3), [-0.60 -0.10], [0 0],'-','Color', [0.35 0.35 0.35],'LineWidth', 0.8);

    %--------------------------------------------------------------
    % Axis labels
    %--------------------------------------------------------------

    xlabel(ax(5), ...
        'Real part \sigma [-]', ...
        'Interpreter', 'tex', ...
        'FontName', 'Arial', ...
        'FontSize', 16, ...
        'FontWeight', 'bold');

    ylabel(ax(3), ...
        'Imaginary part \nu [-]', ...
        'Interpreter', 'tex', ...
        'FontName', 'Arial', ...
        'FontSize', 16, ...
        'FontWeight', 'bold');

    % Text annotations matching the slide exactly
    text(ax(1), -0.52, 2.04, 'High frequency eigenvalue', 'FontSize', 10, 'Color', [0.3 0.3 0.3], 'FontName', 'Arial');
    text(ax(2), -0.52, 1.04, 'Collective and difference eigenvalue', 'FontSize', 10, 'Color', [0.3 0.3 0.3], 'FontName', 'Arial');
    text(ax(3), -0.52, 0.04, 'Low frequency eigenvalue', 'FontSize', 10, 'Color', [0.3 0.3 0.3], 'FontName', 'Arial');
    text(ax(4), -0.52, -1.04, 'Collective and difference eigenvalue', 'FontSize', 10, 'Color', [0.3 0.3 0.3], 'FontName', 'Arial');
    text(ax(5), -0.52, -2.04, 'High frequency eigenvalue', 'FontSize', 10, 'Color', [0.3 0.3 0.3], 'FontName', 'Arial');

    % mu = 0.54 split marker text in middle band
    text(ax(3), 0.33, 0.06, '\mu = 0.54', 'FontSize', 10, 'Color', [0.3 0.3 0.3], 'FontName', 'Arial', 'HorizontalAlignment', 'center');

    %--------------------------------------------------------------
    % Axis-break marks
    %--------------------------------------------------------------

    for b = 1:4

        xBreak = left + width;

        gapTop = bottom(b);

        gapBottom = ...
            bottom(b+1) + height(b+1);

        yBreak = (gapTop + gapBottom) / 2;

        annotation(fig, 'line', [xBreak-0.006 xBreak+0.006],[yBreak+0.010 yBreak],'Color', [0.2 0.2 0.2],'LineWidth', 1.2);
      
        annotation(fig, 'line',[xBreak-0.006 xBreak+0.006],[yBreak yBreak-0.010],'Color', [0.2 0.2 0.2], 'LineWidth', 1.2);

    end

    %--------------------------------------------------------------
    % Direction arrows: increasing advance ratio
    %--------------------------------------------------------------

    k1 = max(1, nMu-10);
    k2 = nMu;

    for j = 1:AnzGl

        c = colours(group(j), :);

        for sgn = [-1 1]

            x1 = real(CharEx(k1, j));
            y1 = sgn * imag(CharEx(k1, j));

            x2 = real(CharEx(k2, j));
            y2 = sgn * imag(CharEx(k2, j));

            dx = x2 - x1;
            dy = y2 - y1;

            if hypot(dx, dy) < 1e-12
                continue;
            end

            for b = 1:5

                insideBand = ...
                    y1 >= bandLimits(b, 1) && ...
                    y1 <= bandLimits(b, 2) && ...
                    y2 >= bandLimits(b, 1) && ...
                    y2 <= bandLimits(b, 2);

                if insideBand

                    quiver(ax(b),x1, y1, dx, dy, 0,'Color', c,'LineWidth', 1.3,'MaxHeadSize', 0.8,'HandleVisibility', 'off');

                end

            end

        end

    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Optional saving
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SaveResults = false;
SaveFigure  = false;

if SaveResults

    save('Workspace_FourBladeRotor.mat', ...
        'Auswahl', 'Blatt', ...
        'gamma', 'd2', 'd3', 'd4', 'ebeta', 'nu0', ...
        'mu_paramVec', 'T', ...
        'CharEx', 'CharExMean', ...
        'CharExRe', 'CharExIm', ...
        'CharMult', 'MonodromieAll', 'cond_A');

end

if SaveFigure

    exportgraphics(gcf, ...
        'CharacteristicExponentTrajectories.png', ...
        'Resolution', 300);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Local functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function A = getSystemMatrix( ...
    psi, mu, gamma, d2, d3, d4, ebeta, nu0, Blatt)

x = zeros(2*Blatt, 1);

[~, A] = SchlagDGL( ...
    psi, x, gamma, d2, d3, d4, ...
    mu, ebeta, nu0, Blatt);

end


function dy = fundamentalMatrixRHS( ...
    psi, y, Afun, mu, n)

Phi = reshape(y, n, n);

A = Afun(psi, mu);

dPhi = A * Phi;

dy = dPhi(:);

end


function V = normalizeColumns(V)

columnNorm = sqrt(sum(abs(V).^2, 1));

V = bsxfun(@rdivide, V, max(columnNorm, eps));

end


function assignment = minimumAssignment(C)
% Exact one-to-one assignment using subset dynamic programming.
% No additional toolbox is required.

n = size(C, 1);

assert(size(C, 2) == n, ...
    'The cost matrix must be square.');

assert(n <= 10, ...
    'This implementation supports up to 10 states.');

nStates = 2^n;

bestCost = inf(nStates, 1);
previousMask = zeros(nStates, 1);
selectedColumn = zeros(nStates, 1);

bestCost(1) = 0;

for mask = 0:nStates-1

    occupied = bitget(mask, 1:n);

    row = sum(occupied) + 1;

    if row > n || ~isfinite(bestCost(mask+1))
        continue;
    end

    availableColumns = find(~occupied);

    for col = availableColumns

        newMask = bitset(mask, col, 1);

        candidateCost = ...
            bestCost(mask+1) + C(row, col);

        if candidateCost < bestCost(newMask+1)

            bestCost(newMask+1) = candidateCost;

            previousMask(newMask+1) = mask;

            selectedColumn(newMask+1) = col;

        end

    end

end

assignment = zeros(n, 1);

mask = nStates-1;

for row = n:-1:1

    assignment(row) = selectedColumn(mask+1);

    mask = previousMask(mask+1);

end

end