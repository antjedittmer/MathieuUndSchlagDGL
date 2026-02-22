clear; clc; close all;

%% === FIGURE FOLDER SETUP ===
fDir = 'figureFolder';
if ~isdir(fDir) %#ok<*ISDIR>
    mkdir(fDir);
end
fDirPeters = fullfile(fDir, 'figureFolderPeters');
if ~isdir(fDirPeters)
    mkdir(fDirPeters);
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


%% === 1. LOAD ARNOLD REFERENCE DATA ===
dDirA = fullfile('dataFolder', 'dataFolder_Arnold_Classic_Symmetric_test');
matA = 'STRUTTscheKarte_D1dot5e-01_SW1dot0e-01_unt0.mat';
filenamStrutt = fullfile(dDirA, matA);
useOldData = useOldData == 1 && exist(filenamStrutt,'file') == 2;
pos0 = get(0,'defaultFigurePosition');

if useOldData == 1
    S_A = load(filenamStrutt);
    nu_A = S_A.CharEx(:,1);
    ReA1 = S_A.CharEx(:,3);
    ReA2 = S_A.CharEx(:,4);
    ImA1 = S_A.CharEx(:,8);
    ImA2 = S_A.CharEx(:,8);  % Trusted Im(s_R2) - Arnold reference

    %% ==== 1.a PLOT CHARACTERISTIC EXPONENT ===
    fig1 = figure('Name', 'Characteristic Exponent from Arnold''s Method', 'Color', 'w');
    fig1.Position = [pos0(1), pos0(2)- 0.15*pos0(4), pos0(3), 1.15*pos0(4)];

    scatter(ReA1, ImA1, 40, nu_A, 'o', 'DisplayName','s_{R1} = \sigma_1 + i \omega_1');
    hold on;
    scatter(ReA2, ImA2, 40, nu_A, '*', 'DisplayName','s_{R2}= \sigma_2 + i \omega_2');
    grid on;

    xlabel('Real part characteristic exponent \sigma')
    ylabel('Pos. imag. part characteristic exponent \omega')

    title('Real vs. Imag. Parts Char. Exponents (Arnolds''s Method)')
    set(findall(gcf,'-property','FontSize'),'FontSize',12)
    cb = colorbar;
    cb.Label.String = 'Amplication factor \nu_c^2';

    legend('Location','southoutside','Orientation','horizontal','FontSize',12);
    legend boxoff

    pngname = sprintf('real_vs_imaginary_exponents_Mathieu_Arnold.png');
    pngfile = fullfile(fDirPeters, pngname);
    saveas(fig1, pngfile);
    % fprintf('Figure saved: %s\n', pngfile);
    % fprintf('Arnold: %d points ✓\n', length(nu_A));

end

%% === PRE-ALLOCATION ===
growth_rate = zeros(length(nu_vals), 1);
frequency_imag = zeros(length(nu_vals), 1);
composite_freq = zeros(length(nu_vals), 1);
participation_data = zeros(length(nu_vals), length(m_range));
branch_freqs_all = zeros(length(nu_vals), length(m_range));
real_all = zeros(length(nu_vals), 2);
total_mag_all = zeros(length(nu_vals), 1);
buffer.Pos = 0;
buffer.Neg = 0;


%% === MAIN COMPUTATION LOOP ===

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


    %idx = 2;
    eta_mode = eta_vals(idx);
    v_mode = V(:, idx);
    eta_prev = eta_mode;
    real_all(k,:) = real(eta_vals(idx2))';

    % Check if we are in a "bubble": Real parts are not identical
    isBubble = abs(real_all(k,1) - real_all(k,2) ) > eps;

    % 'Bubble counting
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
    frequency_imag(k) = mod(Eig.ImagCorrected,0.5); %+ m_bubble(k);


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

%% === PLOT comparison vs amplification factor ===
nc = size(unique(lines,'rows'),1); % Numbers of different colors: nc = 7
lines0 = lines;
cl = [lines0(1:nc,:); 0*ones(1,3); 0.5*ones(1,3)]; % 7 line colors, black and grey

fig = figure('Name', 'Characteristic Exponent: Peters vs Arnold', 'Color', 'w');
fig.Position = [pos0(1), pos0(2)- 0.3*pos0(4), pos0(3), 1.5*pos0(4)];
tiledlayout(4,1,'TileSpacing','tight');

% Subplot 1: Composite Frequency (Peters solid) vs Arnold (dashed) + Growth Rate
nexttile; %subplot(4,1,1);
yyaxis left
if useOldData == 1
    plot(nu_A, ImA2, '-', 'LineWidth', 1.3, 'Color', 'black', 'DisplayName', '\omega (Arnold)');
else
    plot(nu_vals, frequency_imag+ m_bubble,...
        'LineWidth', 1.3, 'Color', 'black', 'DisplayName', '\omega (Arnold)');
end
hold on;

plot(nu_vals, composite_freq, '-.', 'LineWidth', 1.3, 'Color', cl(1,:), 'DisplayName', '\omega (Peters)');
ylabel('\omega = Im(s_R)');
yyaxis right
% plot(nu_A, ReA1, '-', 'LineWidth', 1.3, 'Color', 0.5*ones(1,3) , 'DisplayName', '\sigma (Arnold)');
plot(nu_vals, growth_rate, '-.', 'LineWidth', 1.3, 'Color', cl(2,:), 'DisplayName', '\sigma (Peters)');
ylabel('\sigma = Re(s_R)');
title({['Mathieu ODE: x''''(t)+ 2Dx''(t) +(\nu_0^2 + \nu_c^2 cos(t))x(t) = 0,'...
    sprintf(' D = %2.2f', D),', \nu_0 = \nu_c'] ...
    'Real and Imaginary Part of Characteristic Exponent s_r'});
legend('Location','best');
grid on;

% Subplot 2: Difference (Peters - Arnold)
nexttile; %subplot(4,1,2);
ImA2Interp = frequency_imag+ m_bubble;
if useOldData == 1
    ImA2Interp = (interp1(nu_A, ImA2, nu_vals))';
end
plot(nu_vals, composite_freq - ImA2Interp, 'Color', cl(1,:), 'LineWidth', 1.3, 'DisplayName', 'Peters - Arnold');
ylabel('Δω');
title('Difference Imaginary Part of Characteristic Exponent');
legend('Location','best');
grid on;

% Subplot 3: Frequency Branches
lsCell = {'-','--','-.',':','--'};
nexttile; %subplot(4,1,3);
for idx = 1: size(branch_freqs_all,2)
    idxLs = mod(idx-1,length(lsCell)) + 1;
    idxC = mod(idx-1,length(cl)) + 1;
    %idxLs = floor((idx-1)/length(lsCell)) +1;
    plot(nu_vals, branch_freqs_all(:,idx), 'LineWidth', 1.3, 'color',cl(idxC,:),'LineStyle',lsCell{idxLs});
    hold on;
end
ylabel('Branch Freqs');
title('Individual Frequency Branches');
legend(arrayfun(@(m) sprintf('m=%d', m_range(m)), 1:length(m_range), 'UniformOutput', false), 'Location','northeastoutside');
grid on;

% Subplot 4: Harmonic Participation
nexttile; %subplot(4,1,4);
%set(gca, 'ColorOrder', jet(length(m_range)));
for idx = 1: size(branch_freqs_all,2)
    idxLs = mod(idx-1,length(lsCell)) + 1;
    idxC = mod(idx-1,length(cl)) + 1;
    plot(nu_vals, participation_data(:,idx), 'LineWidth', 1, 'color',cl(idxC,:),'LineStyle',lsCell{idxLs});
    hold on;
end
xlabel('Amplification factor \nu^2_c');
ylabel('Weight');
title('Harmonic Participation');
legend(arrayfun(@(m) sprintf('m=%d', m_range(m)), 1:length(m_range), 'UniformOutput', false), 'Location','northeastoutside');
grid on;

% === SAVE FIGURE ===
pngname = sprintf('compare_char_exp_freq_participation_w1dot0.png');
pngfile = fullfile(fDirPeters, pngname);
saveas(fig, pngfile);
fprintf('Figure saved: %s\n', pngfile);

%% ===

fig2 = figure('Name', 'Characteristic Exponent from Peters'' Method', 'Color', 'w');
fig2.Position = [pos0(1), pos0(2)- 0.15*pos0(4), pos0(3), 1.15*pos0(4)];

scatter(real_all(:,2), composite_freq, 40, nu_vals, 'o', 'DisplayName','\sigma_1 + i |\omega_1|');
hold on;
scatter(real_all(:,1), composite_freq, 40, nu_vals, '*', 'DisplayName',' \sigma_2 + i |\omega_2|');
grid on;

xlabel('Real part characteristic exponent \sigma')
ylabel('Pos. imag. part characteristic exponent \omega')

title('Real vs. Imag. Parts Char. Exponents (Peters'' Method)')
set(findall(gcf,'-property','FontSize'),'FontSize',12)
cb = colorbar;
cb.Label.String = 'Amplication factor \nu_c^2';

lgd = legend('Location','southoutside','Orientation','horizontal','FontSize',12);
legend boxoff

pngname = sprintf('real_vs_imaginary_exponents_Mathieu_Peters.png');
pngfile = fullfile(fDirPeters, pngname);
saveas(fig2, pngfile);


%% ===


fig3 = figure('Name', 'Characteristic Exponent from Arnold''s Method (500 points)', 'Color', 'w');
fig3.Position = [pos0(1), pos0(2)- 0.15*pos0(4), pos0(3), 1.15*pos0(4)];

scatter(real_all(:,2), frequency_imag+ m_bubble, 40, nu_vals, 'o', 'DisplayName','\sigma_1 + i |\omega_1|');
hold on;
scatter(real_all(:,1), frequency_imag+ m_bubble, 40, nu_vals, '*', 'DisplayName','\sigma_2 + i |\omega_2|');
grid on;

xlabel('Real part characteristic exponent \sigma')
ylabel('Pos. imag. part characteristic exponent \omega')

title('Real vs. Imag. Parts Char. Exponents (Arnold''s Method)')
set(findall(gcf,'-property','FontSize'),'FontSize',12)
cb = colorbar;
cb.Label.String = 'Amplication factor \nu_c^2';

lgd = legend('Location','southoutside','Orientation','horizontal','FontSize',12);
legend boxoff

if useOldData == 0
    pngname = sprintf('real_vs_imaginary_exponents_Mathieu_Arnold.png');
    pngfile = fullfile(fDirPeters, pngname);
    saveas(fig3, pngfile);
end


%% 1. Process Modal Participation (Vectorized)

% Define thresholds and constants
gap_threshold = 0.01;
peak_prominence = 0.05;
allModes = 0;

% Sort all rows descending to find dominant branches
[sortedM, idxSortM] = sort(participation_data, 2, 'descend');

% Initialize output matrices
sortedIndex = sortedM(:, 1:2);
%sortedIndexI = sort(idxSortM(:, 1:2), 2); % Standardize ID order

% Mask rows where one mode clearly dominates (the 'else' case)
if allModes == 1
    isNotNearlyEqual = (sortedIndex(:,1) - sortedIndex(:,2)) >= gap_threshold;
    sortedIndex(isNotNearlyEqual, 2)  = 0;
    %sortedIndexI(isNotNearlyEqual, 2) = 0;
end

%% 2. Identify Branch Transitions
% Get peaks of the frequency plot (Peters' composite frequency)
[pksC, locsC] = findpeaks(composite_freq);

% Find transition peaks from participation gap dynamics
% We use the derivative of the difference between top modes to find shifts
diffSorted = [0; diff(sortedIndex(:,2) - sortedIndex(:,1))];
[~, locsD] = findpeaks(diffSorted, 'MinPeakProminence', peak_prominence);

% Detect where bubble counting (m_bubble) changes for comparison
idx_m_change = [false; diff(m_bubble) ~= 0];

%% 3. Generate Winding Number from Modal Participation
m_modpart_raw = zeros(size(m_bubble));
m_modpart_raw(locsD(1)) = 0.5; % Trigger first shift from participation
m_modpart_raw(locsC)    = 0.5; % Trigger subsequent shifts from freq peaks
m_modpart = cumsum(m_modpart_raw);

%% 4. Visualization
pos0 = get(0,'defaultFigurePosition');

cl = lines;
fig4 = figure('Name', 'Floquet Branch Tracking Analysis');
fig4.Position = [pos0(1), pos0(2)- 0.15*pos0(4), pos0(3), 1.4*pos0(4)];


% Subplot 1: Frequency and Peak Detection
subplot(4,1,1)
plot(nu_vals, composite_freq, 'Color', cl(1,:), 'DisplayName', '\omega (Peters)');
hold on;
plot(nu_vals(locsC), pksC, 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'Freq Peaks');
plot(nu_vals(idx_m_change), composite_freq(idx_m_change), 'ro', 'DisplayName', 'm-bubble change');
ylabel('Frequency \omega'); grid on; axis tight;
legend('Location', 'best');

% Subplot 2: Raw Modal Participation
subplot(4,1,2)
plot(nu_vals, participation_data);
hold on;
plot(nu_vals, max(participation_data, [], 2), 'r--', 'LineWidth', 1.3);
plot(nu_vals, sortedIndex(:,2), 'b--', 'LineWidth', 1.3);
ylabel('Mod. Part.'); grid on; axis tight;

% Subplot 3: Participation Gap & Transition Logic
subplot(4,1,3)
plot(nu_vals, sortedIndex(:,1),'r--', 'DisplayName', 'Mode 1');
hold on;
plot(nu_vals, sortedIndex(:,2), 'b--', 'DisplayName', 'Mode 2');
plot(nu_vals, diffSorted, 'k', 'DisplayName', '\Delta (Mode 1 - Mode 2)');
plot(nu_vals(locsD), diffSorted(locsD),  'ko', 'MarkerFaceColor', 'k', 'DisplayName', '\Delta Peaks');
ylabel('Branch Dynamics'); grid on; axis tight;
legend('Location', 'best','NumColumns',2);

% Subplot 4: Winding Number Comparison
subplot(4,1,4)
plot(nu_vals, m_bubble, 'Color', cl(1,:), 'LineWidth', 1.3, 'DisplayName', 'm bubble');
hold on;
plot(nu_vals, m_modpart, '--', 'Color', cl(2,:), 'LineWidth', 1.3, 'DisplayName', 'm Peters');
ylabel('Winding No. m'); grid on; axis tight;
xlabel('Amplification factor \nu');
legend('Location', 'best');

function  [Eig, buffer] = correctImagValues(Eig, buffer)
% correctImagValues korrigiert die Imaginaerwerte fuer stetigen Verlauf
% Inputs
% - Eig: Struct mit Vektor Eig.Imag mit zwei Werten der Imaginaerteile
% - buffer struct: Puffer mit letzten Werten fuer Maximum und Minimum
% Outputs
% - Eig: Struct mit angehaengtem, korrigierten Imaginaerteilen
% - buffer struct: Puffer ueberschrieben mit neuen Werten fuer Maximum und Minimum
% Zwei Checks fuer das korrekte Format
if nargin~= 2
    error('Two inputs are expected: The current imaginary part of the eigenvalues and the buffer with the last values');
end
if max(size(Eig.Imag))~= 2 || min(size(Eig.Imag))~= 1
    error('The current imaginary parts of an eigenvalue pair is expected');
end
% Sortiere Imaginaerteil
Eig.ImagSort = sort(Eig.Imag);
% Imaginaeranteil kontinuierlich steigend oder fallend
tmp = Eig.ImagSort(2);
tmpNeg = Eig.ImagSort(1);

% Initialize buffer.Pos/buffer.Neg if they don't exist, though they should be
% initialized in the main loop to 0
if ~isfield(buffer, 'Pos'), buffer.Pos = 0; end
if ~isfield(buffer, 'Neg'), buffer.Neg = 0; end

if buffer.Pos <= tmp || (abs(tmp) < 10^-5) % Wert uebernehmen
    Eig.ImagCorrected = tmp; % steigender pos. Wert
    Eig.ImagCorrectedNeg = tmpNeg;  % fallender neg. Wert
    buffer.Pos = 0; % Removed based on original logic, these were only placeholders
    buffer.Neg = 0;
else  % 'korrigierten' Wert nehmen fuer kontinuierlichen Verlauf
    % Wird der 'else'-Zweig getriggert, ist der buffer
    % bei 90 Grad (1.57 rad). Der positive Imaginaerteil ist damit
    % die Summe aus 180 Grad und dem negativen Winkel, dessen
    % Wert von -90 Grad zu 0 Grad laeuft.
    Eig.ImagCorrected = 2*buffer.Pos  + tmpNeg; % korrigierter pos. Wert
    Eig.ImagCorrectedNeg  = 2*buffer.Neg + tmp; % korrigierter neg. Wert
end
% Buffer ueberschreiben mit aktuellem Wert
buffer.Pos = max(tmp,buffer.Pos); % Maximum speichern
buffer.Neg = min(tmpNeg,buffer.Neg); % Minimum speichern
% Anmerkung: hier kammen leider andere Werte raus als bei atanh
% Eig.ImagAngle = 1/T * angle(imag(eP)./real(eP));
end
