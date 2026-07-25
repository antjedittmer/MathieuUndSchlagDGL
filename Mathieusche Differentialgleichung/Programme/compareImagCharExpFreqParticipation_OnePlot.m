clear; clc; close all;

%% === FIGURE FOLDER SETUP ===
fDir = 'figureFolder';
if ~isdir(fDir) %#ok<*ISDIR>
    mkdir(fDir);
end
fDirPeters = fullfile(fDir,'figureFolderPeters_SVG');
if ~isdir(fDirPeters)
    mkdir(fDirPeters);
end

dDirMain = fullfile('dataFolder/');
if exist(dDirMain,"dir") ~= 7
    mkdir(dDirMain);
end

dDirPeters = fullfile(dDirMain,'dataFolderOnePlot');
if ~isdir(dDirPeters)
    mkdir(dDirPeters);
end

%% === GLOBAL PARAMETERS ===
Omega = 1;
T = 2*pi/Omega;
nu_vals = linspace(0, 9, 500);
m_range = -4:4;
N_FFT = 2048;
D = 0.15;
x0 = eye(2);
useOldData = 0;

loadData = 1;

%% === 1. LOAD ARNOLD REFERENCE DATA (optional) ===
dDirA = fullfile('dataFolder', 'dataFolder_Arnold_Classic_Symmetric_test');
matA = 'STRUTTscheKarte_D1dot5e-01_SW1dot0e-01_unt0.mat';
filenamStrutt = fullfile(dDirA, matA);
useOldData = useOldData == 1 && exist(filenamStrutt,'file') == 2;
if useOldData == 1
    S_A = load(filenamStrutt);
    nu_A = S_A.CharEx(:,1);
    ImA2 = S_A.CharEx(:,8);  % Trusted Im(s_R2) - Arnold reference
end

%% === PRE-ALLOCATION ===

dataDirPeters = fullfile(dDirPeters,'mathieuFloquetCombined.mat');

if loadData == 1 && exist(dataDirPeters,'file') == 2

    load(dataDirPeters, 'nu_vals', 'D', 'm_range', ...
        'real_all', 'branch_freqs_all', 'participation_data', ...
        'omega_arnold', 'composite_freq', 'freq_argmax', ...
        'sortedIndex', 'diffSorted', 'pksC', 'locsC', 'pksD', 'locsD', ...
        'm_bubble', 'm_modpart', 'nuChange', 'idx_m_change');

else

    growth_rate = zeros(length(nu_vals), 1);
    frequency_imag = zeros(length(nu_vals), 1);
    composite_freq = zeros(length(nu_vals), 1);
    participation_data = zeros(length(nu_vals), length(m_range));
    branch_freqs_all = zeros(length(nu_vals), length(m_range));
    real_all = zeros(length(nu_vals), 2);
    total_mag_all = zeros(length(nu_vals), 1);
    buffer.Pos = 0;
    buffer.Neg = 0;

    %% === MAIN COMPUTATION LOOP (unchanged) ===
    m_bubble = zeros(length(nu_vals), 1);

    for k = 1:length(nu_vals)
        nu = nu_vals(k);

        % Solve for Monodromy Matrix
        ode_mat = @(t, x) [0, 1; -(nu + nu*cos(Omega*t)), -2*D] * reshape(x, 2, 2);
        [~, sol_raw] = ode45(@(t, x) reshape(ode_mat(t, x), 4, 1), [0, T], reshape(x0, 4, 1));
        Phi_T = reshape(sol_raw(end, :), 2, 2);

        % Extract Floquet Exponents (Peters method)
        [V, L_mat] = eig(Phi_T);
        eta_vals = log(diag(L_mat)) / T;

        Eig.Imag = imag(eta_vals);
        [Eig, buffer] = correctImagValues(Eig, buffer);

        [~, idx2] = sort(imag(eta_vals),'descend');
        idx = idx2(1);

        eta_mode = eta_vals(idx);
        v_mode = V(:, idx);
        real_all(k,:) = real(eta_vals(idx2))';

        % Check if we are in a "bubble": Real parts are not identical
        isBubble = abs(real_all(k,1) - real_all(k,2) ) > eps;

        % Bubble counting
        if k == 1 % The first vector value is left at 0
            % Do nothing
        elseif isBubble && ~wasBubble
            % If "bubble" opens, allow m to follow the ODE winding number
            m_bubble(k) = m_bubble(k-1) + 0.5;
        else
            % Otherwise leave m at last value
            m_bubble(k) = m_bubble(k-1);
        end
        wasBubble = isBubble;

        growth_rate(k) = real(eta_mode);
        frequency_imag(k) = mod(Eig.ImagCorrected,0.5);

        % Extract Periodic Part P(t) via FFT for participation
        t_fft = linspace(0, T, N_FFT+1); t_fft(end) = [];
        sol_obj = ode45(@(t, x) reshape(ode_mat(t, x), 4, 1), [0, T], reshape(x0, 4, 1));
        Phi_t_steps = deval(sol_obj, t_fft);

        Q_t = zeros(N_FFT, 1);
        for j = 1:N_FFT
            Phi_curr = reshape(Phi_t_steps(:,j), 2, 2);
            P_t = Phi_curr * v_mode * exp(-eta_mode * t_fft(j));
            Q_t(j) = P_t(1);
        end

        C = fftshift(fft(Q_t)/N_FFT);
        freqs_fft_m = (-(N_FFT/2) : (N_FFT/2-1));

        % Frequency and Participation Calculation
        for i_m = 1:length(m_range)
            m = m_range(i_m);
            [~, f_idx] = min(abs(freqs_fft_m - m));
            participation_data(k, i_m) = abs(C(f_idx));
            branch_freqs_all(k, i_m) = abs(m*Omega + imag(eta_mode));
        end

        % Normalize weights and compute composite frequency
        weights = participation_data(k, :);
        total_mag = sum(weights);
        total_mag_all(k) = total_mag;

        weights = weights / total_mag;
        participation_data(k, :) = weights;
        composite_freq(k) = sum(branch_freqs_all(k, :) .* weights);
    end

    %% === POST-PROCESSING FOR COMBINED FIGURE ===

    % Reference frequency: imaginary part of the char. exponent (Arnold)
    omega_arnold = frequency_imag + m_bubble;
    if useOldData == 1
        omega_arnold = (interp1(nu_A, ImA2, nu_vals))';
    end

    % Frequency of the branch with dominant (argmax) participation
    [sortedM, idxSortM] = sort(participation_data, 2, 'descend');
    iMax = idxSortM(:,1);
    freq_argmax = branch_freqs_all(sub2ind(size(branch_freqs_all), ...
        (1:length(nu_vals))', iMax));

    % Dominant two harmonics ("Mode 1" and "Mode 2") and their gap dynamics
    sortedIndex = sortedM(:, 1:2);
    diffSorted = [0; diff(sortedIndex(:,2) - sortedIndex(:,1))];

    % Peaks of the composite (Peters) frequency
    [pksC, locsC] = findpeaks(composite_freq);

    % Peaks of the participation gap dynamics
    peak_prominence = 0.01;
    [pksD, locsD] = findpeaks(diffSorted, 'MinPeakProminence', peak_prominence);

    % Winding number reconstructed from participation/frequency peaks
    m_modpart_raw = zeros(size(m_bubble));
    if ~isempty(locsD)
        m_modpart_raw(locsD(1)) = 0.5; % Trigger first shift from participation
    end
    m_modpart_raw(locsC) = 0.5;        % Trigger subsequent shifts from freq peaks
    m_modpart = cumsum(m_modpart_raw);

    % m-bubble change locations
    idx_m_change = [false; diff(m_bubble) ~= 0];
    nuChange = nu_vals(idx_m_change);   % positions of the vertical marker lines


    save(dataDirPeters, 'nu_vals', 'D', 'm_range', ...
        'real_all', 'branch_freqs_all', 'participation_data', ...
        'omega_arnold', 'composite_freq', 'freq_argmax', ...
        'sortedIndex', 'diffSorted', 'pksC', 'locsC', 'pksD', 'locsD', ...
        'm_bubble', 'm_modpart', 'nuChange', 'idx_m_change');
end

%% === DOMINANT BRANCH INDEX (three-stage tie-break) ===
% Recomputed here so it also exists when data is loaded from the .mat file.
%
% Inside locked regions the twin branches m and -m-1 are DOUBLY degenerate:
% their participations are (near-)identical AND their frequencies coincide
% exactly (e.g. f(0) = omega = 0.5 and f(-1) = 1 - omega = 0.5). A plain
% argmax and a frequency-based tie-break therefore both pick an arbitrary
% twin. Resolution in three stages, each only active when the previous one
% is genuinely ambiguous:
%   Stage 1: candidates = branches within tolTie of the max participation
%   Stage 2: among them, keep those with (near-)largest branch frequency,
%            f(m) = |m| + sign(m)*omega  (with m = 0 -> omega)
%   Stage 3: among frequency-tied twins, pick the branch whose frequency
%            is closest to omega(Arnold); this fixes the identity inside
%            the lock and hands over to the continuing branch at exit.
% Stage 3 only ever chooses among branches whose participation AND
% frequency agree within tolerance, i.e. where the omega_max VALUE is
% unaffected and only the color label was arbitrary.
tolTie  = 0.01;     % participation tolerance for "nearly identical"
tolFreq = 0.02;     % frequency tolerance for "tied" branch frequencies
n = length(nu_vals);
omega0 = branch_freqs_all(:, m_range == 0);   % |Im(s_R)| principal, in [0, 0.5]

signM = sign(m_range);
signM(m_range == 0) = 1;                      % m = 0 -> frequency omega itself
tieFreq = abs(m_range) + omega0 * signM;      % n x length(m_range)

iMax = zeros(n,1);
for k = 1:n
    % Stage 1: candidates by participation
    vmax = max(participation_data(k,:));
    cand = find(participation_data(k,:) >= vmax - tolTie);
    % Stage 2: among them, keep those with (near-)largest frequency
    fC = tieFreq(k, cand);
    candF = cand(fC >= max(fC) - tolFreq);
    % Stage 3: among frequency-tied twins, pick the branch closest to
    % omega(Arnold)
    [~, jj] = min(abs(tieFreq(k, candF) - omega_arnold(k)));
    iMax(k) = candF(jj);
end

% Recompute the argmax frequency from the tie-broken index so panels 2
% and 4 show the same branch identity (values coincide where branches
% are degenerate, so this only fixes identity, not information)
freq_argmax = branch_freqs_all(sub2ind(size(branch_freqs_all), (1:n)', iMax));

%% === COMBINED FIGURE: 5 stacked axes ===
nc = size(unique(lines,'rows'),1); % Numbers of different colors: nc = 7
lines0 = lines;
cl = [lines0(1:nc,:); 0*ones(1,3); 0.5*ones(1,3)]; % 7 line colors, black, grey
clLight = 1 - 0.60*(1 - cl);  % branch colors: 60 % toward full color (less washed out)
lsCell = {'-','--','-.','-','--'};
clVert = [0 0 0.55];          % dark blue for the m-change vertical lines

pos0 = get(0,'defaultFigurePosition');
fig = figure('Name', 'Floquet Analysis: Combined Diagnostics', 'Color', 'w');
fig.Position = [pos0(1), pos0(2) - 0.45*pos0(4), pos0(3), 1.7*pos0(4)];
tlo = tiledlayout(5,1,'TileSpacing','tight','Padding','compact');
axList = gobjects(5,1);

% --- Axis 1: Real parts of BOTH exponents (bubbles visible) ---
axList(1) = nexttile;
plot(nu_vals, real_all(:,1), '-',  'LineWidth', 1, 'Color', 'b', ...
    'DisplayName', '\sigma_1 = Re(s_{R1})');
hold on;
plot(nu_vals, real_all(:,2), '--', 'LineWidth', 1, 'Color', 'b', ...
    'DisplayName', '\sigma_2 = Re(s_{R2})');
ylabel('\sigma = Re(s_R)');

title({['Mathieu ODE: x''''(t)+ 2Dx''(t) +(\nu_0^2 + \nu_c^2 cos(t))x(t) = 0,',...
    sprintf(' D = %2.2f', D)] ...
    });

legend('Location','eastoutside','Box','off'); grid on; axis tight;
set(gca, 'XTickLabel', []);

% --- Axis 2: Frequency branches + Arnold, centroid, argmax, peaks ---
axList(2) = nexttile;
hBr = gobjects(length(m_range),1);
for idx = 1:size(branch_freqs_all,2)
    idxLs = mod(idx-1,length(lsCell)) + 1;
    idxC  = mod(idx-1,size(cl,1)) + 1;
    hBr(idx) = plot(nu_vals, branch_freqs_all(:,idx), 'LineWidth', 1, ...
        'Color', clLight(idxC,:), 'LineStyle', lsCell{idxLs}, ...
        'DisplayName', sprintf('m = %d', m_range(idx)));
    hold on;
end

hAm = plot(nu_vals, freq_argmax, '.', 'LineWidth', 2 ,'Color', 0.5*ones(3,1), 'MarkerSize',15,...
    'HandleVisibility','off');
hA  = plot(nu_vals, omega_arnold, '-',  'LineWidth', 1.3, 'Color', 'b', ...
    'HandleVisibility','off');
hC  = plot(nu_vals, composite_freq, '--', 'LineWidth', 1.3, 'Color', 'k', ...
    'HandleVisibility','off');
hP  = plot(nu_vals(locsC), pksC, 'ko', 'MarkerFaceColor', 'k', ...
    'MarkerSize', 5, 'HandleVisibility','off');
ylabel('Frequency \omega');
% Legend: only the m-branches; overlay curves are labeled by in-plot text
legend(hBr, 'Location','eastoutside','Box','off'); grid on; axis tight;
set(gca, 'XTickLabel', []);

% In-plot labels for the overlay curves, with tight white backgrounds
bcVec = [1 1 1 0.7];
text(6.35, interp1(nu_vals, omega_arnold, 8.35) + 0.65, ...
    '\omega_1 = Im(s_{R1})', 'Color', 'b', 'FontSize', 11, ...
    'BackgroundColor', bcVec , 'Margin', 0.1);
text(3.30, interp1(nu_vals, freq_argmax, 3.60) - 0.77, ...
    '\omega(\phi_{max})', 'Color', 0.4*ones(1,3), 'FontSize', 13, ...
    'FontWeight', 'bold', 'BackgroundColor', bcVec , 'Margin', 0.1);
text(nu_vals(locsC(3)) - 2, pksC(3) + 0.75, ...
    'freq. peaks \omega_{mean} ', 'Color', 'k', 'FontSize', 11, ...
    'BackgroundColor', bcVec , 'Margin', 0.1);
text(nu_vals(locsC(1)) - 0.6, pksC(1) + 0.75, ...
    '\omega_{mean}', 'Color', 'k', 'FontSize', 11, 'FontWeight', 'bold', ...
    'BackgroundColor', bcVec ,'Margin', 0.1);

% --- Axis 3: Differences w.r.t. Im(s_R) (Arnold) ---
axList(3) = nexttile;

plot(nu_vals, freq_argmax - omega_arnold, '.', 'LineWidth', 1.7, ...
    'Color',  0.5*ones(1,3), 'MarkerSize',15, 'DisplayName', '\omega_{max} - \omega_1');

hold on;

plot(nu_vals, composite_freq - omega_arnold, '--', 'LineWidth', 1.5, ...
    'Color', 'k', 'DisplayName', '\omega_{mean} - \omega_1');

yline(0, 'Color', 0.5*ones(1,3), 'HandleVisibility','off');
ylabel('Difference \Delta\omega');
legend('Location','eastoutside','Box','off'); grid on; axis tight;
set(gca, 'XTickLabel', []);

% --- Axis 4: Harmonic participation, dominant branch highlighted ---
axList(4) = nexttile;
hPart = gobjects(length(m_range),1);
for idx = 1:size(participation_data,2)
    idxLs = mod(idx-1,length(lsCell)) + 1;
    idxC  = mod(idx-1,size(cl,1)) + 1;
    hPart(idx) = plot(nu_vals, participation_data(:,idx), 'LineWidth', 1, ...
        'Color', clLight(idxC,:), 'LineStyle', lsCell{idxLs}, ...
        'DisplayName', sprintf('m = %d', m_range(idx)));
    hold on;
end

% DOMINANT branch highlighted: wherever branch idx has the highest
% participation (three-stage tie-break), its curve segment is drawn solid
% in its own bright branch color -> directly explains the omega_max line
% in panel 2. NaN elsewhere breaks the line into segments automatically.
domCurve = nan(length(nu_vals), length(m_range));
for idx = 1:length(m_range)
    mask = (iMax == idx);
    domCurve(mask, idx) = participation_data(mask, idx);
end
for idx = 1:length(m_range)
    if any(~isnan(domCurve(:,idx)))
        idxC = mod(idx-1,size(cl,1)) + 1;
        plot(nu_vals, domCurve(:,idx), '-', 'LineWidth', 1.7, ...
            'Color', cl(idxC,:), 'HandleVisibility','off');
    end
end

% Gap dynamics between the two dominant harmonics and its peaks
hD  = plot(nu_vals, diffSorted, 'k-.', 'LineWidth', 1, ...
    'HandleVisibility','off');
hDp = plot(nu_vals(locsD), diffSorted(locsD), 'ko', 'MarkerFaceColor', 'k', ...
    'MarkerSize', 5, 'HandleVisibility','off');

ylabel('Mod. Part. \phi');
% Legend: only the m-branches (thin lines); the thick bright segments use
% the same colors, marking where each branch is the dominant one
legend(hPart, 'Location','eastoutside','Box','off'); grid on; axis tight;
set(gca, 'XTickLabel', []);

% In-plot labels (adjust coordinates to taste)
text(0.45, 0.85, 'bold: dominant branch \phi_{max}', 'Color', 'k', ...
    'FontSize', 10, 'BackgroundColor', bcVec, 'Margin', 0.1);
text(1, 0.67, '\Delta = \phi_{max} - \phi_{2nd}', 'Color', 'k', ...
    'FontSize', 10, 'BackgroundColor', bcVec, 'Margin', 1);
if ~isempty(locsD)
    text(nu_vals(locsD(3)) + 0.5, diffSorted(locsD(3)) + 0.35, ...
        'peaks \Delta', 'Color', 'k', 'FontSize', 10, ...
        'BackgroundColor', bcVec, 'Margin', 1);
end

% --- Axis 5: Winding numbers ---
axList(5) = nexttile;
plot(nu_vals, m_bubble, '-', 'Color', 'b', 'LineWidth', 1, ...
    'DisplayName', 'm bubble');
hold on;
plot(nu_vals, m_modpart, '--', 'Color', 'k', 'LineWidth', 1, ...
    'DisplayName', 'm Peters');
ylabel('Add. factor m');
xlabel('Amplification factor \nu_c^2 = \nu_0^2');
legend('Location','eastoutside','Box','off'); grid on; axis tight;

% Explanation of the vertical lines (placed once, in the last panel)
text(2.5, 0.25*max(m_bubble), 'vertical lines: change of m', ...
    'Color', clVert, 'FontSize', 10, 'BackgroundColor', 'w', 'Margin', 1);


% --- Dark blue vertical lines at every m-change, in ALL panels ---
% Drawn as ordinary line objects (not xline/ConstantLine) so that text
% boxes can be stacked ON TOP of them; dotted style
for iAx = 1:5
    yl = ylim(axList(iAx));                 % freeze current limits
    for xv = nuChange
        line(axList(iAx), [xv xv], yl, 'LineStyle', ':', ...
            'Color', clVert, 'LineWidth', 1, 'HandleVisibility', 'off');
    end
    ylim(axList(iAx), yl);                  % restore (lines sit exactly on limits)
    % bring all text boxes of this panel in front of the vertical lines
    % ch = axList(iAx).Children;
    % isTxt = arrayfun(@(h) isa(h, 'matlab.graphics.primitive.Text'), ch);
    % axList(iAx).Children = [ch(isTxt); ch(~isTxt)];
end

% Create arrow
annotation(gcf,'arrow',[0.305714285714286 0.291428571428571],...
    [0.695188340807175 0.673766816143498]);

% Create arrow
annotation(gcf,'arrow',[0.319285714285714 0.292857142857143],...
    [0.311219730941704 0.291479820627803]);

% Create arrow
annotation(gcf,'arrow',[0.234285714285714 0.160714285714286],...
    [0.355381165919283 0.304372197309417]);


% Link x-axes for synchronized zooming
linkaxes(axList, 'x');
xlim(axList(1), [nu_vals(1), nu_vals(end)]);

set(findall(fig, '-property', 'FontSize'), 'FontSize', 11)

%% === SAVE FIGURE ===
basename = 'combined_char_exp_freq_participation_winding';
svgfile = fullfile(fDirPeters, [basename, '.svg']);
pngfile = fullfile(fDirPeters, [basename, '.png']);
print(fig, svgfile, '-dsvg');
print(fig, pngfile, '-dpng', '-r300');
fprintf('Figure saved: %s\n', svgfile);

%% === EXPORT DATA TABLE ===
table4Excel = table;
table4Excel.nu_vals = nu_vals';
table4Excel.char_exp_real1 = real_all(:,1);
table4Excel.char_exp_real2 = real_all(:,2);
table4Excel.char_exp_imag_arnold = omega_arnold;
table4Excel.char_exp_imag_centroid = composite_freq;
table4Excel.char_exp_imag_argmax = freq_argmax;
table4Excel.diff_centroid_arnold = composite_freq - omega_arnold;
table4Excel.diff_argmax_arnold = freq_argmax - omega_arnold;
table4Excel.m_bubble = m_bubble;
table4Excel.m_modpart = m_modpart;

colNames = cell(1, length(m_range));
for i = 1:length(m_range)
    m = m_range(i);
    if m < 0
        colNames{i} = sprintf('freq_branch_m_neg_%d', abs(m));
    else
        colNames{i} = sprintf('freq_branch_m_pos_%d', m);
    end
end
tmpTableFreq = array2table(branch_freqs_all, 'VariableNames', colNames);
colNamesPhi = strrep(colNames,'freq_branch','harm_part');
tmpTableHarmPart = array2table(participation_data, 'VariableNames', colNamesPhi);

table4ExcelAll = [table4Excel, tmpTableFreq, tmpTableHarmPart];
writetable(table4ExcelAll, fullfile(fDirPeters, [basename, '.xlsx']));

%% === LOCAL FUNCTION (unchanged) ===
function  [Eig, buffer] = correctImagValues(Eig, buffer)
% correctImagValues korrigiert die Imaginaerwerte fuer stetigen Verlauf
% Inputs
% - Eig: Struct mit Vektor Eig.Imag mit zwei Werten der Imaginaerteile
% - buffer struct: Puffer mit letzten Werten fuer Maximum und Minimum
% Outputs
% - Eig: Struct mit angehaengtem, korrigierten Imaginaerteilen
% - buffer struct: Puffer ueberschrieben mit neuen Werten fuer Maximum und Minimum
if nargin~= 2
    error('Two inputs are expected: The current imaginary part of the eigenvalues and the buffer with the last values');
end
if max(size(Eig.Imag))~= 2 || min(size(Eig.Imag))~= 1
    error('The current imaginary parts of an eigenvalue pair is expected');
end
% Sortiere Imaginaerteil
Eig.ImagSort = sort(Eig.Imag);
tmp = Eig.ImagSort(2);
tmpNeg = Eig.ImagSort(1);

if ~isfield(buffer, 'Pos'), buffer.Pos = 0; end
if ~isfield(buffer, 'Neg'), buffer.Neg = 0; end

if buffer.Pos <= tmp || (abs(tmp) < 10^-5) % Wert uebernehmen
    Eig.ImagCorrected = tmp; % steigender pos. Wert
    Eig.ImagCorrectedNeg = tmpNeg;  % fallender neg. Wert
    buffer.Pos = 0;
    buffer.Neg = 0;
else  % 'korrigierten' Wert nehmen fuer kontinuierlichen Verlauf
    Eig.ImagCorrected = 2*buffer.Pos  + tmpNeg; % korrigierter pos. Wert
    Eig.ImagCorrectedNeg  = 2*buffer.Neg + tmp; % korrigierter neg. Wert
end
buffer.Pos = max(tmp,buffer.Pos); % Maximum speichern
buffer.Neg = min(tmpNeg,buffer.Neg); % Minimum speichern
end