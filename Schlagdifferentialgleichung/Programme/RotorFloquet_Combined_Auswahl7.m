clc; clear variables; close all;

%% === FOLDER SETUP ===
excelDir = 'excelDir';
if ~isfolder(excelDir), mkdir(excelDir); end
figDir = 'figDir';
if ~isfolder(figDir), mkdir(figDir); end

%% === ROTOR SELECTION ===
% This combined-figure method (correctImagValues with Pos/Neg buffers,
% "bubble" winding number, dominant-branch participation) is written for
% a 2-state characteristic-exponent pair, matching Auswahl 7
% (Einzelblattkoordinaten, Blatt = 1, AnzGl = 2) directly.
Auswahl = 7;
Blatt   = 1;
rotorDescription = 'Einzelblattkoordinaten im rotierenden System';
AnzGl   = 2;
konstant = 0;

fprintf('Running: Auswahl %d - %s\n', Auswahl, rotorDescription);

%% === PARAMETERS ===
Parameter = readtable('Parameter.xlsx', 'Range', 'C4:I29', 'ReadVariableNames', false);
Par = table2array(Parameter);

ebeta = Par(8,  Auswahl);
gamma = Par(13, Auswahl);
d2    = Par(17, Auswahl);
d3    = Par(18, Auswahl);
d4    = Par(19, Auswahl);
nu0   = Par(20, Auswahl);

T  = 2*pi;
SW = 0.1;
MuMin = 0; MuMax = 10;
mu_vals = MuMin:SW:MuMax;
n = numel(mu_vals);

m_range = -4:4;
mAll    = numel(m_range);
N_FFT   = 4096;
x0      = eye(AnzGl);
flapStateIndex = 2;   % Blatt == 1 case

options2 = odeset('RelTol', 1e-10, 'AbsTol', 1e-12, 'MaxStep', 1e-3);
t_fft = linspace(0, T, N_FFT+1); t_fft(end) = [];
freqs_fft_m = (-(N_FFT/2) : (N_FFT/2-1));

%% === PRE-ALLOCATION ===
growth_rate       = zeros(n, 1);
frequency_imag    = zeros(n, 1);
composite_freq    = zeros(n, 1);
participation_data = zeros(n, mAll);
branch_freqs_all  = zeros(n, mAll);
real_all          = zeros(n, 2);
cond_A            = nan(n, 1);
m_bubble          = zeros(n, 1);

buffer.Pos = 0;
buffer.Neg = 0;
wasBubble  = false;

%% === MAIN COMPUTATION LOOP ===
for k = 1:n
    mu_param = mu_vals(k);

    sol_cell = cell(1, AnzGl);
    Phi_T    = zeros(AnzGl);
    for j = 1:AnzGl
        if konstant == 1
            sol = ode45(@(psi, x) SchlagDGLkonstant(psi, x, gamma, d2, d3, d4, ...
                mu_param, ebeta, nu0, Blatt), [0, T], x0(:, j), options2);
        else
            sol = ode45(@(psi, x) SchlagDGL(psi, x, gamma, d2, d3, d4, ...
                mu_param, ebeta, nu0, Blatt), [0, T], x0(:, j), options2);
        end
        sol_cell{j} = sol;
        Phi_T(:, j) = deval(sol, T);
    end
    cond_A(k) = cond(Phi_T);

    % --- Extract Floquet Exponents ---
    [V, L_mat] = eig(Phi_T);
    eta_vals = log(diag(L_mat)) / T;

    Eig.Imag = imag(eta_vals);
    [Eig, buffer] = correctImagValues(Eig, buffer);

    [~, idx2] = sort(imag(eta_vals), 'descend');
    idx = idx2(1);

    eta_mode = eta_vals(idx);
    v_mode   = V(:, idx);
    real_all(k, :) = real(eta_vals(idx2))';

    % "Bubble" = the pair has split into two distinct real growth rates
    % instead of remaining a complex-conjugate (oscillatory) pair.
    isBubble = abs(real_all(k, 1) - real_all(k, 2)) > eps;

    if k == 1
        % m_bubble(1) stays 0
    elseif isBubble && ~wasBubble
        m_bubble(k) = m_bubble(k-1) + 0.5;
    else
        m_bubble(k) = m_bubble(k-1);
    end
    wasBubble = isBubble;

    growth_rate(k)    = real(eta_mode);
    frequency_imag(k) = mod(Eig.ImagCorrected, 0.5);

    % --- Harmonic participation via FFT of the periodic part ---
    Q_t = complex(zeros(N_FFT, 1));
    for it = 1:N_FFT
        Phi_t = zeros(AnzGl);
        for j = 1:AnzGl
            Phi_t(:, j) = deval(sol_cell{j}, t_fft(it));
        end
        P_t     = Phi_t * v_mode * exp(-eta_mode * t_fft(it));
        Q_t(it) = P_t(flapStateIndex);
    end

    C = fftshift(fft(Q_t) / N_FFT);

    for i_m = 1:mAll
        m = m_range(i_m);
        [~, f_idx] = min(abs(freqs_fft_m - m));
        participation_data(k, i_m) = abs(C(f_idx));
        branch_freqs_all(k, i_m)   = abs(m + imag(eta_mode));
    end

    weights   = participation_data(k, :);
    total_mag = sum(weights);
    if total_mag > 1e-12
        weights = weights / total_mag;
    else
        weights = zeros(size(weights));
    end
    participation_data(k, :) = weights;
    composite_freq(k)        = sum(branch_freqs_all(k, :) .* weights);
end

%% === POST-PROCESSING ===
omega_track = frequency_imag + m_bubble;

[sortedM, idxSortM] = sort(participation_data, 2, 'descend');
sortedIndex = sortedM(:, 1:2);
diffSorted  = [0; diff(sortedIndex(:,2) - sortedIndex(:,1))];

[pksC, locsC] = findpeaks(composite_freq);

peak_prominence = 0.01;
[pksD, locsD] = findpeaks(diffSorted, 'MinPeakProminence', peak_prominence);

m_modpart_raw = zeros(size(m_bubble));
if ~isempty(locsD)
    m_modpart_raw(locsD(1)) = 0.5;
end
m_modpart_raw(locsC) = 0.5;
m_modpart = cumsum(m_modpart_raw);

idx_m_change = [false; diff(m_bubble) ~= 0];
muChange     = mu_vals(idx_m_change);

% Reliability guard - see condReliableMax note below
condReliableMax = 1e9;
firstUnreliableMu = NaN;
if any(cond_A > condReliableMax)
    firstUnreliableMu = mu_vals(find(cond_A > condReliableMax, 1, 'first'));
    fprintf(['Warning: cond(Phi) exceeds %.0e starting at mu = %.2f - ' ...
             'results beyond this mu are floating-point noise, not physics.\n'], ...
             condReliableMax, firstUnreliableMu);
end

%% === DOMINANT BRANCH INDEX (argmax with tie-break among locked twins) ===
distPart = [(mAll-1)/2 - 1 + 0.5 : -1 : 0.5, 0:(mAll-1)/2];  % height per branch, m_range order
if numel(distPart) ~= mAll
    % Fallback for even-length m_range (shouldn't trigger for m_range = -6:6)
    distPart = 0:(mAll-1);
end
tolTie = 0.001;

iMax = zeros(n, 1);
for k = 1:n
    [~, i0] = max(participation_data(k, :));
    cand = find(abs(participation_data(k, :) - participation_data(k, i0)) < tolTie);
    [~, jj] = max(distPart(cand));
    iMax(k) = cand(jj);
end
iMax = iMax(:);

freq_argmax = branch_freqs_all(sub2ind(size(branch_freqs_all), (1:n)', iMax));
m_argmax    = distPart(iMax);
m_argmax    = m_argmax(:);

%% === COMBINED FIGURE: 5 stacked axes ===
nc = size(unique(lines, 'rows'), 1);
lines0 = lines;
cl = [lines0(1:nc, :); 0*ones(1,3); 0.5*ones(1,3)];
clLight = 1 - 0.60*(1 - cl);
lsCell = {'-','--','-.','-','--'};
clVert = [0 0 0.55];

pos0 = get(0, 'defaultFigurePosition');
fig = figure('Name', 'Rotor Floquet Analysis: Combined Diagnostics', 'Color', 'w');
fig.Position = [pos0(1), pos0(2) - 0.45*pos0(4), pos0(3), 1.7*pos0(4)];
tlo = tiledlayout(5, 1, 'TileSpacing', 'tight', 'Padding', 'compact');
axList = gobjects(5, 1);

% --- Axis 1: Real parts of both exponents ---
axList(1) = nexttile;
plot(mu_vals, real_all(:,1), '-',  'LineWidth', 1, 'Color', 'b', ...
    'DisplayName', '\sigma_1 = Re(s_{R1})');
hold on;
plot(mu_vals, real_all(:,2), '--', 'LineWidth', 1, 'Color', 'b', ...
    'DisplayName', '\sigma_2 = Re(s_{R2})');
ylabel('\sigma = Re(s_R)');
title(sprintf('Auswahl %d; Blatt %d; %s', Auswahl, Blatt, rotorDescription));
legend('Location', 'eastoutside', 'Box', 'off'); grid on; axis tight;
set(gca, 'XTickLabel', []);

% --- Axis 2: Frequency branches + tracked frequency, centroid, argmax, peaks ---
axList(2) = nexttile;
hBr = gobjects(mAll, 1);
for idxB = 1:mAll
    idxLs = mod(idxB-1, numel(lsCell)) + 1;
    idxC  = mod(idxB-1, size(cl,1)) + 1;
    hBr(idxB) = plot(mu_vals, branch_freqs_all(:, idxB), 'LineWidth', 1, ...
        'Color', clLight(idxC,:), 'LineStyle', lsCell{idxLs}, ...
        'DisplayName', sprintf('m = %d', m_range(idxB)));
    hold on;
end
plot(mu_vals, freq_argmax, '.', 'Color', 0.5*ones(3,1), 'MarkerSize', 15, 'HandleVisibility', 'off');
plot(mu_vals, omega_track, '-',  'LineWidth', 1.3, 'Color', 'b', 'HandleVisibility', 'off');
plot(mu_vals, composite_freq, '--', 'LineWidth', 1.3, 'Color', 'k', 'HandleVisibility', 'off');
if ~isempty(locsC)
    plot(mu_vals(locsC), pksC, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5, 'HandleVisibility', 'off');
end
ylabel('Frequency \omega');
legend(hBr, 'Location', 'eastoutside', 'Box', 'off'); grid on; axis tight;
set(gca, 'XTickLabel', []);

% --- Axis 3: Differences w.r.t. tracked frequency ---
axList(3) = nexttile;
plot(mu_vals, freq_argmax - omega_track, '.', 'LineWidth', 1.7, ...
    'Color', 0.5*ones(1,3), 'MarkerSize', 15, 'DisplayName', '\omega_{max} - \omega_{track}');
hold on;
plot(mu_vals, composite_freq - omega_track, '--', 'LineWidth', 1.5, ...
    'Color', 'k', 'DisplayName', '\omega_{mean} - \omega_{track}');
yline(0, 'Color', 0.5*ones(1,3), 'HandleVisibility', 'off');
ylabel('Difference \Delta\omega');
legend('Location', 'eastoutside', 'Box', 'off'); grid on; axis tight;
set(gca, 'XTickLabel', []);

% --- Axis 4: Harmonic participation, dominant branch highlighted ---
axList(4) = nexttile;
hPart = gobjects(mAll, 1);
for idxB = 1:mAll
    idxLs = mod(idxB-1, numel(lsCell)) + 1;
    idxC  = mod(idxB-1, size(cl,1)) + 1;
    hPart(idxB) = plot(mu_vals, participation_data(:, idxB), 'LineWidth', 1, ...
        'Color', clLight(idxC,:), 'LineStyle', lsCell{idxLs}, ...
        'DisplayName', sprintf('m = %d', m_range(idxB)));
    hold on;
end

domCurve = nan(n, mAll);
for idxB = 1:mAll
    mask = (iMax == idxB);
    maskExt = mask | [false; mask(1:end-1)] | [mask(2:end); false];
    domCurve(maskExt, idxB) = participation_data(maskExt, idxB);
end
for idxB = 1:mAll
    if any(~isnan(domCurve(:, idxB)))
        idxC = mod(idxB-1, size(cl,1)) + 1;
        plot(mu_vals, domCurve(:, idxB), '-', 'LineWidth', 1.7, ...
            'Color', cl(idxC,:), 'HandleVisibility', 'off');
    end
end

plot(mu_vals, diffSorted, 'k-.', 'LineWidth', 1, 'HandleVisibility', 'off');
if ~isempty(locsD)
    plot(mu_vals(locsD), diffSorted(locsD), 'ko', 'MarkerFaceColor', 'k', ...
        'MarkerSize', 5, 'HandleVisibility', 'off');
end
ylabel('Mod. Part. \phi');
legend(hPart, 'Location', 'eastoutside', 'Box', 'off'); grid on; axis tight;
set(gca, 'XTickLabel', []);

% --- Axis 5: Winding numbers ---
axList(5) = nexttile;
plot(mu_vals, m_bubble, '-', 'Color', 'b', 'LineWidth', 1, 'DisplayName', 'm bubble');
hold on;
plot(mu_vals, m_modpart, '--', 'Color', 'k', 'LineWidth', 1, 'DisplayName', 'm Peters');

mArgCurve = nan(n, mAll);
for idxB = 1:mAll
    mArgCurve(iMax == idxB, idxB) = distPart(idxB);
end
for idxB = 1:mAll
    if any(~isnan(mArgCurve(:, idxB)))
        idxC = mod(idxB-1, size(cl,1)) + 1;
        plot(mu_vals, mArgCurve(:, idxB), '.', 'Color', cl(idxC,:), 'MarkerSize', 12, ...
            'HandleVisibility', 'off');
    end
end
plot(NaN, NaN, '.', 'Color', 0.4*ones(1,3), 'MarkerSize', 12, 'DisplayName', 'm(\phi_{max})');

ylabel('Add. factor m');
xlabel('Advance ratio \mu (-)');
legend('Location', 'eastoutside', 'Box', 'off'); grid on; axis tight;

% --- Vertical lines at every m-bubble change, in all panels ---
for iAx = 1:5
    yl = ylim(axList(iAx));
    for xv = muChange
        line(axList(iAx), [xv xv], yl, 'LineStyle', ':', ...
            'Color', clVert, 'LineWidth', 1, 'HandleVisibility', 'off');
    end
    ylim(axList(iAx), yl);
end

% --- Mark the numerically-unreliable region (cond(Phi) too large) ---
if ~isnan(firstUnreliableMu)
    for iAx = 1:5
        xline(axList(iAx), firstUnreliableMu, 'r:', 'LineWidth', 1.2, ...
            'HandleVisibility', 'off');
    end
end

linkaxes(axList, 'x');
xlim(axList(1), [mu_vals(1), mu_vals(end)]);
set(findall(fig, '-property', 'FontSize'), 'FontSize', 11);

%% === SAVE FIGURE ===
basename = sprintf('RotorFloquet_Combined_Auswahl%d', Auswahl);
svgfile = fullfile(figDir, [basename, '.svg']);
pngfile = fullfile(figDir, [basename, '.png']);
print(fig, svgfile, '-dsvg');
print(fig, pngfile, '-dpng', '-r300');
fprintf('Figure saved: %s\n', svgfile);

%% === EXPORT DATA TABLE ===
table4Excel = table;
table4Excel.mu               = mu_vals';
table4Excel.char_exp_real1   = real_all(:,1);
table4Excel.char_exp_real2   = real_all(:,2);
table4Excel.char_exp_imag_track = omega_track;
table4Excel.char_exp_imag_centroid = composite_freq;
table4Excel.char_exp_imag_argmax   = freq_argmax;
table4Excel.diff_centroid_track = composite_freq - omega_track;
table4Excel.diff_argmax_track   = freq_argmax - omega_track;
table4Excel.m_bubble  = m_bubble;
table4Excel.m_modpart = m_modpart;
table4Excel.m_argmax  = m_argmax;
table4Excel.m_sel_index = m_range(iMax)';
table4Excel.cond_Phi   = cond_A;

colNames = cell(1, mAll);
for i = 1:mAll
    m = m_range(i);
    if m < 0
        colNames{i} = sprintf('freq_branch_m_neg_%d', abs(m));
    else
        colNames{i} = sprintf('freq_branch_m_pos_%d', m);
    end
end
tmpTableFreq = array2table(branch_freqs_all, 'VariableNames', colNames);
colNamesPhi  = strrep(colNames, 'freq_branch', 'harm_part');
tmpTableHarmPart = array2table(participation_data, 'VariableNames', colNamesPhi);

table4ExcelAll = [table4Excel, tmpTableFreq, tmpTableHarmPart];
writetable(table4ExcelAll, fullfile(excelDir, [basename, '.xlsx']));

%% === LOCAL FUNCTION (unchanged from reference) ===
function [Eig, buffer] = correctImagValues(Eig, buffer)
% correctImagValues korrigiert die Imaginaerwerte fuer stetigen Verlauf
if nargin ~= 2
    error('Two inputs are expected: the current imaginary parts and the buffer.');
end
if max(size(Eig.Imag)) ~= 2 || min(size(Eig.Imag)) ~= 1
    error('The current imaginary parts of an eigenvalue pair are expected.');
end

Eig.ImagSort = sort(Eig.Imag);
tmp    = Eig.ImagSort(2);
tmpNeg = Eig.ImagSort(1);

if ~isfield(buffer, 'Pos'), buffer.Pos = 0; end
if ~isfield(buffer, 'Neg'), buffer.Neg = 0; end

if buffer.Pos <= tmp || (abs(tmp) < 1e-5)
    Eig.ImagCorrected    = tmp;
    Eig.ImagCorrectedNeg = tmpNeg;
    buffer.Pos = 0;
    buffer.Neg = 0;
else
    Eig.ImagCorrected    = 2*buffer.Pos + tmpNeg;
    Eig.ImagCorrectedNeg = 2*buffer.Neg + tmp;
end
buffer.Pos = max(tmp, buffer.Pos);
buffer.Neg = min(tmpNeg, buffer.Neg);
end