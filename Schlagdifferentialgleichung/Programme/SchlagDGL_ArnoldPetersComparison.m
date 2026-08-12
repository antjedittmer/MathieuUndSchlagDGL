%% HELICOPTER ROTOR FLAPWISE STABILITY & HARMONIC PARTICIPATION ANALYSIS
% Analyzes characteristic exponents via Floquet theory (Monodromy matrix)
% and maps flapwise harmonic branch tracking across advance ratios.

clc; clear; close all;

%% =========================================================================
%% 1. INITIALIZATION & CONFIGURATION
%% =========================================================================
SW       = 0.1;
t0       = 0.0;
T        = 2 * pi;

excelDir = 'excelDir';
figDir   = 'figDir';

if ~isfolder(excelDir)
    mkdir(excelDir);
end
if ~isfolder(figDir)
    mkdir(figDir);
end

N_FFT             = 4096;
m_range           = -6:6;
n_harmonics       = length(m_range);
minPeakProminence = 0.02;
minPeakDistanceMu = 0.4;
plotswitch        = 1;

%% =========================================================================
%% 2. ROTOR SELECTION CONFIGURATION
%% =========================================================================
AuswahlInfo = {
    1, 3, '3-Blatt-Rotor, see-saw';
    2, 3, '3-Blatt-Rotor, voll gelenkig';
    3, 4, '4-Blatt-Rotor, voll gelenkig';
    4, 5, '5-Blatt-Rotor, voll gelenkig';
    5, 3, '3-Blatt-Rotor, gelenk-/lagerlos';
    6, 4, '4-Blatt-Rotor, gelenk-/lagerlos';
    7, 1, 'Einzelblattkoordinaten im rotierenden System'
    };

Auswahl          = 7;
Blatt            = AuswahlInfo{Auswahl, 2};
rotorDescription = AuswahlInfo{Auswahl, 3};

AnzGl  = 2 * Blatt;
b1     = AnzGl;
nPairs = Blatt;

fprintf('Running: Auswahl %d — %s\n', Auswahl, rotorDescription);

%% =========================================================================
%% 3. PARAMETER LOADING (EXCEL ENTRY)
%% =========================================================================
Parameter = readtable('Parameter.xlsx', 'Range', 'C4:I29');
Par       = table2array(Parameter);

rho   = Par(26, 1); %#ok<NASGU>
ebeta = Par(8, Auswahl);
gamma = Par(13, Auswahl);
d2    = Par(17, Auswahl);
d3    = Par(18, Auswahl);
d4    = Par(19, Auswahl);
nu0   = Par(20, Auswahl);

MuMin  = 0;
MuMax  = 10;
mu_vec = MuMin:SW:MuMax;
nMu    = numel(mu_vec);

%% =========================================================================
%% 4. FLOQUET CHARACTERISTIC EXPONENTS SYSTEM PASS
%% =========================================================================
Diagonal    = eye(AnzGl);
Monodromie  = zeros(AnzGl);
CharMult    = zeros(nMu, AnzGl + 1);
CharExRe    = zeros(nMu, AnzGl);
CharExIm    = zeros(nMu, AnzGl);
CharExImRaw = zeros(nMu, AnzGl);
CharEx      = zeros(nMu, AnzGl*6);
cond_A      = nan(nMu, 1);

buffer.Pos  = zeros(nPairs, 1);

idx                 = 1;
nShift              = 1;
prevAbsDiffCharMult = 0;

options2 = odeset('RelTol', 1e-10, 'AbsTol', 1e-12, 'MaxStep', 1e-3);

for mu_param = mu_vec

    for k = 1:AnzGl
        sol = ode45(@(psi, x) SchlagDGL(psi, x, gamma, d2, d3, d4, ...
            mu_param, ebeta, nu0, Blatt), [t0, T], Diagonal(:, k), options2);
        Monodromie(:, k) = deval(sol, T);
    end

    cond_A(idx) = cond(Monodromie);

    charMult = eig(Monodromie).';
    [~, idxSort] = sort(real(charMult));
    charMultSort = charMult(idxSort);

    CharMult(idx, :) = [charMultSort, mu_param];
    if mu_param == 0
        Im0 = 1/T * angle(charMultSort);
    end

   % Charakteristische Exponenten
            Eig.Real = 1/T * log(abs(charMultSort));
            Eig.Imag = 1/T * atan(imag(charMultSort) ./ real(charMultSort));

            % correctImagValues paarweise
            ImagCorrected    = zeros(1, nPairs);
            ImagCorrectedNeg = zeros(1, nPairs);

    for idxC = 1:nPairs
        idxVec         = 2*idxC-1 : 2*idxC;
        EigTemp.Real   = Eig.Real(idxVec);
        EigTemp.Imag   = Eig.Imag(idxVec);
        bufferTemp.Pos = buffer.Pos(idxC);

        [EigTemp, bufferTemp] = correctImagValuesLocal(EigTemp, bufferTemp);

        Eig.Real(idxVec)     = EigTemp.Real;
        Eig.Imag(idxVec)     = EigTemp.Imag;
        Eig.ImagSort(idxVec) = EigTemp.ImagSort;
        ImagCorrected(idxC)    = EigTemp.ImagCorrected;
        ImagCorrectedNeg(idxC) = EigTemp.ImagCorrectedNeg;
        buffer.Pos(idxC)     = bufferTemp.Pos;
    end
    Eig.ImagCorrected    = ImagCorrected;
    Eig.ImagCorrectedNeg = ImagCorrectedNeg;

    absDiffCharMult = abs(CharMult(idx, 1)) > abs(CharMult(idx, b1));

    if idx > 1 && abs(prevAbsDiffCharMult - absDiffCharMult) > 0.1
        nShift = nShift + 0.5;
       
    end

    prevAbsDiffCharMult = absDiffCharMult;

    nAdd    = nShift * 2 * pi / T;
    nAddVec = [-nAdd * ones(1, nPairs), nAdd * ones(1, nPairs)];

    ImagEigSortN = [Eig.ImagCorrectedNeg, Eig.ImagCorrected] + nAddVec;

            CharEx(idx,:)     = [Eig.Real, Eig.Imag, Eig.ImagSort, ImagEigSortN, ...
                                  real(charMultSort), imag(charMultSort)];
            CharExRe(idx,:)    = Eig.Real;
            CharExIm(idx,:)    = Eig.Imag + nAddVec;
            CharExImRaw(idx,:) = Eig.Imag;

    idx = idx + 1;
end

%% =========================================================================
%% 5. DATA POST-PROCESSING & BRANCH-ALIGNED TRACKING
%% =========================================================================
CharExRe1 = CharExRe(:, 1);
CharExRe2 = CharExRe(:, Blatt + 1);

CharExRePos      = max(CharExRe1, CharExRe2);
CharExRe1_NegIdx = abs(CharExRePos - CharExRe1) > eps;

offset      = CharExRe1(1);
CharExReNeg = -(CharExRePos - offset) + offset;

CharExRe1Cor = CharExRe1;
CharExRe1Cor(CharExRe1_NegIdx) = CharExReNeg(CharExRe1_NegIdx);

CharExRe2Cor = CharExRe2;
CharExRe2Cor(~CharExRe1_NegIdx) = CharExReNeg(~CharExRe1_NegIdx);

CharExIm1 = CharExIm(:, 1);
CharExIm2 = CharExIm(:, Blatt + 1);

% Use tolerance, not exact equality
tolMask   = 1e-10;
isBranch1 = abs(CharExRe1Cor - CharExRe1) < tolMask;

% Raw tracked imaginary branches aligned with the same branch logic
CharExIm1RawTracked = zeros(nMu,1);
CharExIm2RawTracked = zeros(nMu,1);

CharExIm1RawTracked(isBranch1)  = CharExIm1(isBranch1);
CharExIm2RawTracked(isBranch1)  = -CharExIm1(isBranch1);

CharExIm1RawTracked(~isBranch1) = -CharExIm2(~isBranch1);
CharExIm2RawTracked(~isBranch1) =  CharExIm2(~isBranch1);

% Corrected branches use the same aligned representation
CharExIm1Cor = CharExIm1RawTracked;
CharExIm2Cor = CharExIm2RawTracked;

%% =========================================================================
%% 6. FLOQUET EXCURSION DIAGRAM PLOTTING
%% =========================================================================
cl   = lines;
fig1 = figure(Auswahl); clf;
pos0 = get(groot, 'DefaultFigurePosition');

set(fig1, 'Position', [pos0(1:2), 1.1 * pos0(3), pos0(4)])

tl = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('Auswahl: %d; Blatt: %d; %s', Auswahl, Blatt, rotorDescription));

% Real Component Subplot
ax2(1) = nexttile;
plot(mu_vec, CharExRe1Cor, '-', 'Color', cl(1,:), 'LineWidth', 1.5); hold on;
plot(mu_vec, CharExRe2Cor, '-', 'Color', cl(2,:), 'LineWidth', 1.5);
plot(mu_vec, CharExRe1, 'k--', 'LineWidth', 1);
plot(mu_vec, CharExRe2, 'k-.', 'LineWidth', 1);
grid on;
ylabel('Real part (-)');
legend('Re(Exp1) cor', sprintf('Re(Exp%d) cor', Blatt+1), ...
    'Re(Exp1) raw', sprintf('Re(Exp%d) raw', Blatt+1), ...
    'Location', 'northeastoutside');

% Imaginary Component Subplot
ax2(2) = nexttile;
plot(mu_vec, CharExIm1Cor, '-', 'Color', cl(1,:), 'LineWidth', 1.5); hold on;
plot(mu_vec, CharExIm2Cor, '-', 'Color', cl(2,:), 'LineWidth', 1.5);
plot(mu_vec, CharExIm1RawTracked, 'k--', 'LineWidth', 1);
plot(mu_vec, CharExIm2RawTracked, 'k-.', 'LineWidth', 1);
grid on;
xlabel('Advance ratio \mu (-)');
ylabel('Imaginary part (-)');
legend('Im(Exp1) cor', sprintf('Im(Exp%d) cor', Blatt+1), ...
    'Im(Exp1) raw tracked', sprintf('Im(Exp%d) raw tracked', Blatt+1), ...
    'Location', 'northeastoutside');

linkaxes(ax2, 'x');

saveas(fig1, fullfile(figDir, sprintf('FloquetExcursion_Auswahl%d.png', Auswahl)));

%% =========================================================================
%% 7. FLAPWISE HARMONIC PARTICIPATION & SPECTRUM ANALYSIS
%% =========================================================================
if Blatt == 1
    flapStateIndex = 2;
else
    flapStateIndex = Blatt + 1;
end

participation_data = zeros(nMu, n_harmonics);
branch_freqs_all   = zeros(nMu, n_harmonics);
composite_freq     = zeros(nMu, 1);
growth_rate        = zeros(nMu, 1);

eta_store        = nan(nMu, 1);
lambda_store     = nan(nMu, 1);
mode_index_store = nan(nMu, 1);

eta_mode_prev = [];

psi_fft      = linspace(0, T, N_FFT + 1);
psi_fft(end) = [];
freq_indices = -N_FFT/2 : N_FFT/2 - 1;

for kMu = 1:nMu
    mu_param = mu_vec(kMu);
    X0       = eye(AnzGl);
    sol_ode  = cell(1, AnzGl);

    for j = 1:AnzGl
        sol_ode{j} = ode45(@(psi, x) SchlagDGL(psi, x, gamma, d2, d3, d4, ...
            mu_param, ebeta, nu0, Blatt), [t0, T], X0(:, j), options2);
    end

    Phi_T = zeros(AnzGl, AnzGl);
    for j = 1:AnzGl
        Phi_T(:, j) = deval(sol_ode{j}, T);
    end

    [V, L]  = eig(Phi_T);
    lambda  = diag(L);
    eta_all = log(lambda) / T;

    if kMu == 1
        [~, mode_idx] = max(imag(eta_all));
    else
        [~, mode_idx] = min(abs(eta_all - eta_mode_prev));
    end

    eta_mode = eta_all(mode_idx);
    v_mode   = V(:, mode_idx);

    [~, idxMax] = max(abs(v_mode));
    v_mode = v_mode / abs(v_mode(idxMax));

    eta_mode_prev = eta_mode;

    eta_store(kMu)        = eta_mode;
    lambda_store(kMu)     = lambda(mode_idx);
    mode_index_store(kMu) = mode_idx;

    growth_rate(kMu) = real(eta_mode);

    Q_t = complex(zeros(N_FFT, 1));
    for it = 1:N_FFT
        Phi_t = zeros(AnzGl, AnzGl);
        for j = 1:AnzGl
            Phi_t(:, j) = deval(sol_ode{j}, psi_fft(it));
        end

        A_t     = Phi_t * v_mode * exp(-eta_mode * psi_fft(it));
        Q_t(it) = A_t(flapStateIndex);
    end

    C_coeffs = fftshift(fft(Q_t) / N_FFT);

    harmonic_magnitudes_raw = zeros(1, n_harmonics);
    for i_m = 1:n_harmonics
        m = m_range(i_m);
        [~, idxFreq] = min(abs(freq_indices - m));
        harmonic_magnitudes_raw(i_m) = abs(C_coeffs(idxFreq));
        branch_freqs_all(kMu, i_m)   = abs(m + imag(eta_mode));
    end

    total_sum = sum(harmonic_magnitudes_raw);
    if total_sum > 1e-12
        participation_data(kMu, :) = harmonic_magnitudes_raw / total_sum;
    end

    weights             = participation_data(kMu, :);
    composite_freq(kMu) = sum(branch_freqs_all(kMu, :) .* weights);
end

%% =========================================================================
%% 8. MODE SORTING & CRITICAL PEAK INTERACTION DETECTION
%% =========================================================================
[part_sorted, mode_rank_idx] = sort(participation_data, 2, 'descend');

mode1_part = part_sorted(:, 1);
mode2_part = part_sorted(:, 2);

mode1_idx = mode_rank_idx(:, 1);
mode2_idx = mode_rank_idx(:, 2);

mode1_m = m_range(mode1_idx).';
mode2_m = m_range(mode2_idx).';

diff12 = [0; abs(diff(mode1_part - mode2_part))];

[pksD, locsD_mu] = findpeaks(diff12, mu_vec, ...
    'MinPeakProminence', minPeakProminence, ...
    'MinPeakDistance', minPeakDistanceMu);

locsD = zeros(size(locsD_mu));
for ii = 1:length(locsD_mu)
    [~, locsD(ii)] = min(abs(mu_vec - locsD_mu(ii)));
end

[pksC, locsC_mu] = findpeaks(composite_freq, mu_vec, ...
    'MinPeakProminence', minPeakProminence, ...
    'MinPeakDistance', minPeakDistanceMu);

locsC = zeros(size(locsC_mu));
for ii = 1:length(locsC_mu)
    [~, locsC(ii)] = min(abs(mu_vec - locsC_mu(ii)));
end

m_factor_raw = zeros(size(mu_vec));
if ~isempty(locsD)
    m_factor_raw(locsD) = m_factor_raw(locsD) + 0.5;
end
if ~isempty(locsC)
    m_factor_raw(locsC) = m_factor_raw(locsC) + 0.5;
end

n        = cumsum(m_factor_raw);
m_peters = n;

%% =========================================================================
%% 9. METRIC PARTICIPATION VISUALIZATION
%% =========================================================================
if plotswitch == 1
    lines0 = lines;
    nc     = size(lines0, 1);
    cl2    = [lines0(1:nc, :); 0 0 0; 0.5 0.5 0.5];
    lsCell = {'-', '--', '-.', ':', '-'};

    fig4 = figure('Name', 'Flapwise Branch Tracking Analysis', 'Color', 'w');
    fig4.Position = [pos0(1), pos0(2) - 0.15 * pos0(4), pos0(3), 1.45 * pos0(4)];
    tiledlayout(4, 1, 'TileSpacing', 'tight', 'Padding', 'compact');

    nexttile;
    plot(mu_vec, composite_freq, 'Color', cl2(1,:), 'LineWidth', 1.3, 'DisplayName', '\omega (Peters)');
    hold on;
    plot(locsC_mu, pksC, 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'Freq Peaks');
    ylabel('Frequency');
    grid on;
    axis tight;
    legend('Location', 'best');
    title(sprintf('Flapwise participation tracking; Auswahl %d; state %d', Auswahl, flapStateIndex));

    nexttile;
    for idxP = 1:size(participation_data, 2)
        idxLs = mod(idxP - 1, length(lsCell)) + 1;
        idxC  = mod(idxP - 1, size(cl2, 1)) + 1;
        plot(mu_vec, participation_data(:, idxP), 'LineWidth', 1.0, ...
            'Color', cl2(idxC, :), 'LineStyle', lsCell{idxLs});
        hold on;
    end
    plot(mu_vec, mode1_part, 'r--', 'LineWidth', 1.3, 'DisplayName', '1st participation');
    plot(mu_vec, mode2_part, 'b--', 'LineWidth', 1.3, 'DisplayName', '2nd participation');
    ylabel('Mod. Part.');
    grid on;
    axis tight;
    title('Raw Harmonic Participation');

    nexttile;
    plot(mu_vec, mode1_part, 'r--', 'DisplayName', 'Mode 1');
    hold on;
    plot(mu_vec, mode2_part, 'b--', 'DisplayName', 'Mode 2');
    plot(mu_vec, diff12, 'k', 'LineWidth', 1.2, 'DisplayName', '\Delta (Mode 1 - Mode 2)');
    plot(locsD_mu, pksD, 'ko', 'MarkerFaceColor', 'k', 'DisplayName', '\Delta Peaks');
    ylabel('Branch Dynamics');
    grid on;
    axis tight;
    legend('Location', 'best', 'NumColumns', 2);
    title('Participation Gap and Peak Detection');

    nexttile;
    plot(mu_vec, n, 'k', 'LineWidth', 1.5, 'DisplayName', 'n');
    hold on;
    plot(mu_vec, m_peters, '--', 'Color', cl2(2,:), 'LineWidth', 1.3, 'DisplayName', 'm Peters');
    ylabel('Factor');
    xlabel('Advance ratio \mu');
    grid on;
    axis tight;
    legend('Location', 'best');
    title('Peak-Triggered Multiplication Factor');

    saveas(fig4, fullfile(figDir, sprintf('FlapwiseBranchTracking_Auswahl%d.png', Auswahl)));
end




% %% =========================================================================
% %% 10. SCATTER PLOT: REAL VS IMAGINARY PARTS ON SAME AXES
% %% =========================================================================
% figScatter = figure('Name', 'Characteristic Exponents in Complex Plane', 'Color', 'w');
% figScatter.Position = [pos0(1), pos0(2)-0.15*pos0(4), pos0(3), 1.15*pos0(4)];
% 
% scatter(CharExRe1Cor(:), CharExIm1Cor(:), 40, mu_vec(:), 'o', ...
%     'DisplayName', 's_{R1} = \sigma_1 + i\omega_1');
% hold on;
% scatter(CharExRe2Cor(:), CharExIm2Cor(:), 40, mu_vec(:), '*', ...
%     'DisplayName', 's_{R2} = \sigma_2 + i\omega_2');
% 
% grid on;
% xlabel('Real part characteristic exponent \sigma')
% ylabel('Pos. imag. part characteristic exponent \omega')
% title(sprintf('Real vs. Imag. Parts Char. Exponents (%s)', rotorDescription))
% set(findall(gcf,'-property','FontSize'),'FontSize',12)
% 
% cb = colorbar;
% cb.Label.String = 'Advance ratio \mu';
% 
% legend('Location','southoutside','Orientation','horizontal','FontSize',12);
% legend boxoff
% 
% saveas(figScatter, fullfile(figDir, sprintf('real_vs_imaginary_exponents_rotor_Auswahl%d.png', Auswahl)));

%% =========================================================================
%% 10. SCATTER PLOTS: PRINCIPAL BRANCH VS TRACKED BRANCH
%% =========================================================================

% Principal-branch exponents for the selected pair
% These are mathematically consistent but not physically unique in frequency.
CharExRe1_principal = CharExRe(:, 1);
CharExRe2_principal = CharExRe(:, Blatt + 1);

% IMPORTANT:
% CharExImRaw contains the principal imaginary part from log(lambda)/T
CharExIm1_principal = CharExImRaw(:, 1);
CharExIm2_principal = CharExImRaw(:, Blatt + 1);

figScatter = figure('Name', 'Characteristic Exponents in Complex Plane', 'Color', 'w');
figScatter.Position = [pos0(1), pos0(2)-0.15*pos0(4), 1.25*pos0(3), 1.05*pos0(4)];

tlSc = tiledlayout(1,2, 'TileSpacing', 'compact', 'Padding', 'compact');

% -------------------------------------------------------------------------
% LEFT: principal branch
% -------------------------------------------------------------------------
axS1 = nexttile;
scatter(CharExRe1_principal(:), CharExIm1_principal(:), 38, mu_vec(:), 'o', ...
    'DisplayName', 'Exp1 principal');
hold on;
scatter(CharExRe2_principal(:), CharExIm2_principal(:), 38, mu_vec(:), '^', ...
    'DisplayName', sprintf('Exp%d principal', Blatt+1));
grid on;
xlabel('Real part characteristic exponent \sigma');
ylabel('Imaginary part \omega (principal branch)');
title('Principal-branch Floquet exponents');
legend('Location','southoutside','Orientation','horizontal');
cb1 = colorbar;
cb1.Label.String = 'Advance ratio \mu';

% -------------------------------------------------------------------------
% RIGHT: tracked branch
% -------------------------------------------------------------------------
axS2 = nexttile;
scatter(CharExRe1Cor(:), CharExIm1Cor(:), 38, mu_vec(:), 'o', ...
    'DisplayName', 'Exp1 tracked');
hold on;
scatter(CharExRe2Cor(:), CharExIm2Cor(:), 38, mu_vec(:), '^', ...
    'DisplayName', sprintf('Exp%d tracked', Blatt+1));
grid on;
xlabel('Real part characteristic exponent \sigma');
ylabel('Imaginary part \omega (tracked branch)');
title('Tracked Floquet branch visualization');
legend('Location','southoutside','Orientation','horizontal');
cb2 = colorbar;
cb2.Label.String = 'Advance ratio \mu';

set(findall(figScatter,'-property','FontSize'),'FontSize',12)

saveas(figScatter, fullfile(figDir, ...
    sprintf('real_vs_imaginary_exponents_rotor_Auswahl%d.png', Auswahl)));
%% =========================================================================
%% 11. FILE EXPORT SYSTEM (.XLSX & .MAT STORAGE)
%% =========================================================================
tableMain = table(mu_vec(:), ...
    CharExRe1Cor(:), CharExRe2Cor(:), CharExIm1Cor(:), CharExIm2Cor(:), ...
    CharExIm1RawTracked(:), CharExIm2RawTracked(:), ...
    composite_freq(:), growth_rate(:), ...
    mode1_part(:), mode2_part(:), diff12(:), ...
    mode1_idx(:), mode2_idx(:), mode1_m(:), mode2_m(:), ...
    n(:), m_peters(:), ...
    'VariableNames', {'mu','CharExRe1Cor','CharExRe2Cor','CharExIm1Cor','CharExIm2Cor', ...
    'CharExIm1RawTracked','CharExIm2RawTracked', ...
    'composite_freq','growth_rate', ...
    'mode1_part','mode2_part','diff12', ...
    'mode1_idx','mode2_idx','mode1_m','mode2_m', ...
    'n','m_peters'});

tableBranch = array2table(branch_freqs_all, ...
    'VariableNames', matlab.lang.makeValidName(strcat("freq_m_", string(m_range))));

tablePart = array2table(participation_data, ...
    'VariableNames', matlab.lang.makeValidName(strcat("part_m_", string(m_range))));

peakTableFreq = table(locsC_mu(:), pksC(:), ...
    'VariableNames', {'mu_peak_freq','peak_freq'});

peakTableDiff = table(locsD_mu(:), pksD(:), ...
    'VariableNames', {'mu_peak_diff','peak_diff'});

finalTable = [tableMain, tableBranch, tablePart];

writetable(finalTable, fullfile(excelDir, sprintf('FlapwiseBranchTracking_Auswahl%d.xlsx', Auswahl)));
writetable(peakTableFreq, fullfile(excelDir, sprintf('FlapwiseFreqPeaks_Auswahl%d.xlsx', Auswahl)));
writetable(peakTableDiff, fullfile(excelDir, sprintf('FlapwiseDiffPeaks_Auswahl%d.xlsx', Auswahl)));

save(fullfile(excelDir, sprintf('FlapwiseBranchTracking_Auswahl%d.mat', Auswahl)), ...
    'finalTable','peakTableFreq','peakTableDiff', ...
    'participation_data','branch_freqs_all','composite_freq','growth_rate', ...
    'part_sorted','mode_rank_idx', ...
    'mode1_part','mode2_part','diff12', ...
    'mode1_idx','mode2_idx','mode1_m','mode2_m', ...
    'n','m_peters', ...
    'locsC_mu','pksC','locsD_mu','pksD', ...
    'mu_vec','m_range','CharExRe1Cor','CharExRe2Cor', ...
    'CharExIm1Cor','CharExIm2Cor','CharExIm1RawTracked','CharExIm2RawTracked', ...
    'Auswahl','Blatt','flapStateIndex');

%% =========================================================================
%% LOCAL FUNCTIONS BLOCK
%% =========================================================================
function [Eig, buffer] = correctImagValuesLocal(Eig, buffer)

if numel(Eig.Imag) ~= 2
    error('Expected exactly one eigenvalue pair.');
end

Eig.ImagSort = sort(Eig.Imag);
tmp    = Eig.ImagSort(2);
tmpNeg = Eig.ImagSort(1);

if buffer.Pos <= tmp || abs(tmp) < 1e-5
    Eig.ImagCorrected    = tmp;
    Eig.ImagCorrectedNeg = tmpNeg;
    buffer.Pos           = 0;
else
    Eig.ImagCorrected    = 2 * buffer.Pos + tmpNeg;
    Eig.ImagCorrectedNeg = tmp;
end

buffer.Pos = max(tmp, buffer.Pos);

end