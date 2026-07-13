clc; clear; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SW       = 0.1;
t0       = 0.0;
T        = 2*pi;
excelDir = 'excelDir';
if ~isfolder(excelDir), mkdir(excelDir); end
figDir = 'figDir';
if ~isfolder(figDir), mkdir(figDir); end

N_FFT             = 4096;
m_range           = -6:6;
minPeakProminence = 0.02;
minPeakDistanceMu = 0.4;
plotswitch        = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Rotor selection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AuswahlInfo = {
    1, 3, '3-Blatt-Rotor, see-saw';
    2, 3, '3-Blatt-Rotor, voll gelenkig';
    3, 4, '4-Blatt-Rotor, voll gelenkig';
    4, 5, '5-Blatt-Rotor, voll gelenkig';
    5, 3, '3-Blatt-Rotor, gelenk-/lagerlos';
    6, 4, '4-Blatt-Rotor, gelenk-/lagerlos';
    7, 1, 'Einzelblattkoordinaten im rotierenden System'
};

Auswahl = 7;
Blatt   = AuswahlInfo{Auswahl,2};
rotorDescription = AuswahlInfo{Auswahl,3};
AnzGl   = 2*Blatt;
b1      = AnzGl;
nPairs  = Blatt;

fprintf('Running: Auswahl %d — %s\n', Auswahl, rotorDescription);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Parameter = readtable('Parameter.xlsx', 'Range', 'C4:I29');
Par       = table2array(Parameter);

rho   = Par(26,1); %#ok<NASGU>
ebeta = Par(8, Auswahl);
gamma = Par(13, Auswahl);
d2    = Par(17, Auswahl);
d3    = Par(18, Auswahl);
d4    = Par(19, Auswahl);
nu0   = Par(20, Auswahl);

MuMin = 0;
MuMax = 10;
mu_vec = MuMin:SW:MuMax;
nMu    = numel(mu_vec);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Characteristic exponents
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Diagonal    = eye(AnzGl);
Monodromie  = zeros(AnzGl);
CharMult    = zeros(nMu, AnzGl+1);
CharExRe    = zeros(nMu, AnzGl);
CharExIm    = zeros(nMu, AnzGl);
CharExImRaw = zeros(nMu, AnzGl);
cond_A      = nan(nMu, 1);
buffer.Pos  = zeros(nPairs, 1);

idx = 1;
nShift = 1;
prevAbsDiffCharMult = 0;

options2 = odeset('RelTol',1e-10,'AbsTol',1e-12,'MaxStep',1e-3);

for mu_param = mu_vec
    for k = 1:AnzGl
        sol = ode45(@(psi,x) SchlagDGL(psi, x, gamma, d2, d3, d4, ...
            mu_param, ebeta, nu0, Blatt), [t0,T], Diagonal(:,k), options2);
        Monodromie(:,k) = deval(sol, T);
    end

    cond_A(idx) = cond(Monodromie);
    charMult    = eig(Monodromie).';
    [~,idxSort] = sort(real(charMult));
    charMultSort = charMult(idxSort);

    CharMult(idx,:) = [charMultSort, mu_param];

    Eig.Real = real(log(charMultSort))/T;
    Eig.Imag = imag(log(charMultSort))/T;

    ImagCorrected    = zeros(1, nPairs);
    ImagCorrectedNeg = zeros(1, nPairs);
    Eig.ImagSort     = zeros(1, AnzGl);

    for idxC = 1:nPairs
        idxVec = 2*idxC-1 : 2*idxC;
        EigTemp.Real   = Eig.Real(idxVec);
        EigTemp.Imag   = Eig.Imag(idxVec);
        bufferTemp.Pos = buffer.Pos(idxC);

        [EigTemp, bufferTemp] = correctImagValuesLocal(EigTemp, bufferTemp);

        Eig.Real(idxVec)     = EigTemp.Real;
        Eig.Imag(idxVec)     = EigTemp.Imag;
        Eig.ImagSort(idxVec) = EigTemp.ImagSort;
        ImagCorrected(idxC)    = EigTemp.ImagCorrected;
        ImagCorrectedNeg(idxC) = EigTemp.ImagCorrectedNeg;
        buffer.Pos(idxC) = bufferTemp.Pos;
    end

    absDiffCharMult = abs(CharMult(idx,1)) > abs(CharMult(idx,b1));
    if idx > 1 && abs(prevAbsDiffCharMult - absDiffCharMult) > 0.1
        nShift = nShift + 0.5;
    end
    prevAbsDiffCharMult = absDiffCharMult;

    nAdd    = nShift * 2*pi / T;
    nAddVec = [-nAdd*ones(1,nPairs), nAdd*ones(1,nPairs)];

    CharExRe(idx,:)    = Eig.Real;
    CharExIm(idx,:)    = Eig.Imag + nAddVec;
    CharExImRaw(idx,:) = Eig.Imag;

    idx = idx + 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Corrected characteristic exponent parts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CharExRe1 = CharExRe(:,1);
CharExRe2 = CharExRe(:,Blatt+1);

CharExRePos      = max(CharExRe1, CharExRe2);
CharExRe1_NegIdx = abs(CharExRePos - CharExRe1) > eps;
offset           = CharExRe1(1);
CharExReNeg      = -(CharExRePos - offset) + offset;

CharExRe1Cor = CharExRe1;
CharExRe1Cor(CharExRe1_NegIdx) = CharExReNeg(CharExRe1_NegIdx);

CharExRe2Cor = CharExRe2;
CharExRe2Cor(~CharExRe1_NegIdx) = CharExReNeg(~CharExRe1_NegIdx);

CharExIm1 = CharExIm(:,1);
CharExIm2 = CharExIm(:,Blatt+1);

CharExIm1Cor = CharExIm1;
CharExIm2Cor = CharExIm2;

for k = 1:nMu
    if CharExRe1Cor(k) == CharExRe1(k)
        CharExIm1Cor(k) =  CharExIm1(k);
        CharExIm2Cor(k) = -CharExIm1(k);
    else
        CharExIm1Cor(k) = -CharExIm2(k);
        CharExIm2Cor(k) =  CharExIm2(k);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Characteristic exponent plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cl = lines;
figure(Auswahl); clf;
pos0 = get(0,"defaultFigurePosition");
set(gcf, "Position", [pos0(1:2), 1.1*pos0(3), pos0(4)])
tl = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
title(tl, sprintf('Auswahl: %d; Blatt: %d; %s', Auswahl, Blatt, rotorDescription));

ax2(1) = nexttile;
plot(mu_vec, CharExRe1Cor,   '-',   'Color', cl(1,:), 'LineWidth', 1.5); hold on;
plot(mu_vec, CharExRe2Cor,   '-',   'Color', cl(2,:), 'LineWidth', 1.5);
plot(mu_vec, CharExRe(:,1),  'k--', 'LineWidth', 1);
plot(mu_vec, CharExRe(:,b1), 'k-.', 'LineWidth', 1);
grid on; ylabel('Real part (-)');
legend('Re(Exp1) cor', sprintf('Re(Exp%d) cor', Blatt), ...
       'Re(Exp1) raw', sprintf('Re(Exp%d) raw', Blatt), ...
       'Location','northeastoutside');

ax2(2) = nexttile;
plot(mu_vec, CharExIm1Cor,   '-',   'Color', cl(1,:), 'LineWidth', 1.5); hold on;
plot(mu_vec, CharExIm2Cor,   '-',   'Color', cl(2,:), 'LineWidth', 1.5);
plot(mu_vec, CharExIm(:,1),  'k--', 'LineWidth', 1);
plot(mu_vec, CharExIm(:,b1), 'k-.', 'LineWidth', 1);
grid on;
xlabel('Advance ratio \mu (-)');
ylabel('Imaginary part (-)');
legend('Im(Exp1) cor', sprintf('Im(Exp%d) cor', Blatt), ...
       'Im(Exp1) raw', sprintf('Im(Exp%d) raw', Blatt), ...
       'Location','northeastoutside');
linkaxes(ax2,'x');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Flapwise harmonic participation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Blatt == 1
    flapStateIndex = 2;
else
    flapStateIndex = Blatt + 1;
end

participation_data = zeros(nMu, length(m_range));
branch_freqs_all   = zeros(nMu, length(m_range));
composite_freq     = zeros(nMu, 1);
growth_rate        = zeros(nMu, 1);
eta_store          = nan(nMu,1);
lambda_store       = nan(nMu,1);
mode_index_store   = nan(nMu,1);

eta_mode_prev = [];

for kMu = 1:nMu
    mu_param = mu_vec(kMu);

    X0 = eye(AnzGl);
    sol_ode = cell(1, AnzGl);

    for j = 1:AnzGl
        sol_ode{j} = ode45(@(psi,x) SchlagDGL(psi, x, gamma, d2, d3, d4, ...
            mu_param, ebeta, nu0, Blatt), [t0, T], X0(:,j), options2);
    end

    Phi_T = zeros(AnzGl, AnzGl);
    for j = 1:AnzGl
        Phi_T(:,j) = deval(sol_ode{j}, T);
    end

    [V, L] = eig(Phi_T);
    lambda = diag(L);
    eta_all = log(lambda) / T;

    if kMu == 1
        [~, mode_idx] = max(imag(eta_all));
    else
        [~, mode_idx] = min(abs(eta_all - eta_mode_prev));
    end

    eta_mode = eta_all(mode_idx);
    v_mode   = V(:, mode_idx);

    [~, idxMax] = max(abs(v_mode));
    v_mode = v_mode / v_mode(idxMax);
    v_mode = v_mode * exp(-1i*angle(v_mode(idxMax)));

    eta_mode_prev = eta_mode;
    eta_store(kMu)      = eta_mode;
    lambda_store(kMu)   = lambda(mode_idx);
    mode_index_store(kMu) = mode_idx;
    growth_rate(kMu)    = real(eta_mode);

    psi_fft = linspace(0, T, N_FFT+1);
    psi_fft(end) = [];

    Q_t = complex(zeros(N_FFT,1));

    for it = 1:N_FFT
        Phi_t = zeros(AnzGl, AnzGl);
        for j = 1:AnzGl
            Phi_t(:,j) = deval(sol_ode{j}, psi_fft(it));
        end
        A_t = Phi_t * v_mode * exp(-eta_mode * psi_fft(it));
        Q_t(it) = A_t(flapStateIndex);
    end

    C_coeffs = fftshift(fft(Q_t) / N_FFT);
    freq_indices = -N_FFT/2 : N_FFT/2 - 1;

    harmonic_magnitudes_raw = zeros(1, length(m_range));
    for i_m = 1:length(m_range)
        m = m_range(i_m);
        [~, idxFreq] = min(abs(freq_indices - m));
        harmonic_magnitudes_raw(i_m) = abs(C_coeffs(idxFreq));
        branch_freqs_all(kMu, i_m)   = abs(m + imag(eta_mode));
    end

    total_sum = sum(harmonic_magnitudes_raw);
    if total_sum > 1e-12
        participation_data(kMu, :) = harmonic_magnitudes_raw / total_sum;
    end

    weights = participation_data(kMu, :);
    composite_freq(kMu) = sum(branch_freqs_all(kMu, :) .* weights);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sort participations and peak detection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[part_sorted, mode_rank_idx] = sort(participation_data, 2, 'descend');
mode1_part = part_sorted(:,1);
mode2_part = part_sorted(:,2);
mode1_idx  = mode_rank_idx(:,1);
mode2_idx  = mode_rank_idx(:,2);

mode1_m = m_range(mode1_idx).';
mode2_m = m_range(mode2_idx).';

diff12 = [0; abs(diff(mode1_part - mode2_part))];

[pksD, locsD_mu] = findpeaks(diff12, mu_vec, ...
    'MinPeakProminence', minPeakProminence, ...
    'MinPeakDistance',   minPeakDistanceMu);

locsD = zeros(size(locsD_mu));
for ii = 1:length(locsD_mu)
    [~, locsD(ii)] = min(abs(mu_vec - locsD_mu(ii)));
end

% n = cumsum(ismember(1:nMu, locsD)).';

[pksC, locsC_mu] = findpeaks(composite_freq, mu_vec, ...
    'MinPeakProminence', minPeakProminence, ...
    'MinPeakDistance',   minPeakDistanceMu);

locsC = zeros(size(locsC_mu));
for ii = 1:length(locsC_mu)
    [~, locsC(ii)] = min(abs(mu_vec - locsC_mu(ii)));
end

m_peters_raw = zeros(size(mu_vec(:)));
if ~isempty(locsD)
    m_peters_raw(locsD(1)) = 0.5;
end
if ~isempty(locsC)
    m_peters_raw(locsC) = m_peters_raw(locsC) + 0.5;
end
m_peters = cumsum(m_peters_raw);

m_modpart_raw = zeros(size(mu_vec));
m_modpart_raw(locsD) = 0.5; % Trigger first shift from participation
m_modpart_raw(locsC) = 0.5; % Trigger subsequent shifts from freq peaks
n =  cumsum(m_modpart_raw);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot participation analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plotswitch == 1
    pos0 = get(0,'defaultFigurePosition');
    lines0 = lines;
    nc = size(lines0,1);
    cl2 = [lines0(1:nc,:); 0 0 0; 0.5 0.5 0.5];
    lsCell = {'-','--','-.',':','-'};

    fig4 = figure('Name', 'Flapwise Branch Tracking Analysis', 'Color', 'w');
    fig4.Position = [pos0(1), pos0(2)-0.15*pos0(4), pos0(3), 1.45*pos0(4)];
    tiledlayout(4,1,'TileSpacing','tight','Padding','compact');

    nexttile
    plot(mu_vec, composite_freq, 'Color', cl2(1,:), 'LineWidth', 1.3, ...
        'DisplayName', '\omega (Peters)');
    hold on;
    plot(locsC_mu, pksC, 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'Freq Peaks');
    ylabel('Frequency');
    grid on; axis tight;
    legend('Location','best');
    title(sprintf('Flapwise participation tracking; Auswahl %d; state %d', ...
        Auswahl, flapStateIndex));

    nexttile
    for idxP = 1:size(participation_data,2)
        idxLs = mod(idxP-1,length(lsCell)) + 1;
        idxC  = mod(idxP-1,size(cl2,1)) + 1;
        plot(mu_vec, participation_data(:,idxP), 'LineWidth', 1.0, ...
            'Color', cl2(idxC,:), 'LineStyle', lsCell{idxLs});
        hold on;
    end
    plot(mu_vec, mode1_part, 'r--', 'LineWidth', 1.3, 'DisplayName', '1st participation');
    plot(mu_vec, mode2_part, 'b--', 'LineWidth', 1.3, 'DisplayName', '2nd participation');
    ylabel('Mod. Part.');
    grid on; axis tight;
    title('Raw Harmonic Participation');

    nexttile
    plot(mu_vec, mode1_part, 'r--', 'DisplayName', 'Mode 1'); hold on;
    plot(mu_vec, mode2_part, 'b--', 'DisplayName', 'Mode 2');
    plot(mu_vec, diff12, 'k', 'LineWidth', 1.2, 'DisplayName', '\Delta (Mode 1 - Mode 2)');
    plot(locsD_mu, pksD, 'ko', 'MarkerFaceColor', 'k', 'DisplayName', '\Delta Peaks');
    ylabel('Branch Dynamics');
    grid on; axis tight;
    legend('Location','best','NumColumns',2);
    title('Participation Gap and Peak Detection');

    nexttile
    plot(mu_vec, n, 'k', 'LineWidth', 1.5, 'DisplayName', 'n'); hold on;
    plot(mu_vec, m_peters, '--', 'Color', cl2(2,:), 'LineWidth', 1.3, ...
        'DisplayName', 'm Peters');
    ylabel('Factor');
    xlabel('Advance ratio \mu');
    grid on; axis tight;
    legend('Location','best');
    title('Peak-Triggered Multiplication Factor');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Export data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tableMain = table(mu_vec(:), ...
    CharExRe1Cor(:), CharExRe2Cor(:), CharExIm1Cor(:), CharExIm2Cor(:), ...
    composite_freq(:), growth_rate(:), ...
    mode1_part(:), mode2_part(:), diff12(:), ...
    mode1_idx(:), mode2_idx(:), mode1_m(:), mode2_m(:), ...
    n(:), m_peters(:), ...
    'VariableNames', {'mu','CharExRe1Cor','CharExRe2Cor','CharExIm1Cor','CharExIm2Cor', ...
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

finalTable = [tableMain tableBranch tablePart];

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
    'mu_vec','m_range','CharExRe1Cor','CharExRe2Cor','CharExIm1Cor','CharExIm2Cor', ...
    'Auswahl','Blatt','flapStateIndex');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Local functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Eig, buffer] = correctImagValuesLocal(Eig, buffer)
    if numel(Eig.Imag) ~= 2
        error('Expected one eigenvalue pair.');
    end

    Eig.ImagSort = sort(Eig.Imag);
    tmp    = Eig.ImagSort(2);
    tmpNeg = Eig.ImagSort(1);

    if buffer.Pos <= tmp || abs(tmp) < 1e-5
        Eig.ImagCorrected    = tmp;
        Eig.ImagCorrectedNeg = tmpNeg;
        buffer.Pos = 0;
    else
        Eig.ImagCorrected    = 2*buffer.Pos + tmpNeg;
        Eig.ImagCorrectedNeg = tmp;
    end

    buffer.Pos = max(tmp, buffer.Pos);
end