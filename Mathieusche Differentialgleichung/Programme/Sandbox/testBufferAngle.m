clc; clear; close all;

%% Lade die Daten
matFolder = fullfile(fileparts(pwd), 'dataFolder', 'dataFolderSTRUTT'); 
matName = fullfile(matFolder, 'STRUTTscheKarte_D1dot5e-01_SW1dot0e-01_unt0.mat'); 

if ~exist(matName, 'file')
    error('File not found: %s', matName);
end

load(matName, 'Char*', 'plotwert*', 'nu*');

% Note: Ensure we use the correct eigenvalues for the branch we want to track
tableCharEx = array2table(CharEx,'VariableNames',...
    {'nu_02','nu_C2','RealEig1','RealEig2','ImagEig1_0','ImagEig2_0',...
    'ImagEig1','ImagEig2','ePidx1','ePidx2','eP1','eP2'});


eP1 = tableCharEx.ePidx1; % eP1 column
eP2 = tableCharEx.ePidx2; % eP1 column
T = 2*pi;

%% 2. Geometric Calculation (Angle + Unwrap)
% Get the total accumulated rotation
rawPhase = angle(eP1);
unwrappedPhase = unwrap(rawPhase);
ImagTotal = abs(unwrappedPhase) / T; % Total monotonic progression

% Calculate the m-factor accumulated at each resonance crossing.
% Using a small epsilon (1e-5) prevents numerical noise from 
% triggering a premature jump at exactly 0.5.
m_factor = floor(ImagTotal * 2 + 1e-5) / 2; 
ImagGeometricSawtooth = ImagTotal - m_factor;

%% 3. Manual Buffer Logic (for comparison)
tmp = CharEx(:, 6); 
tmpNeg = CharEx(:, 5); 
ImagEigManual = nan(size(tmp));
buffer = tmp(1);
ImagEigManual(1) = tmp(1);

for idx = 2:length(tmp)
    if buffer <= tmp(idx) || (abs(tmp(idx)) < 1e-5)
        ImagEigManual(idx) = tmp(idx);
        buffer = 0;
    else
        % Jump logic to match the reset behavior
        ImagEigManual(idx) = 2*buffer + tmpNeg(idx);
    end
    buffer = max(tmp(idx), buffer);
end

%% 4. Results & Plotting
figure(12); clf;
set(gcf, 'Name', 'Corrected Geometric Sawtooth');

subplot(2,1,1)
plot(ImagEigManual, 'LineWidth', 2, 'DisplayName', 'Manual Buffer Logic'); hold on;
plot(ImagGeometricSawtooth, '--', 'LineWidth', 1.5, 'DisplayName', 'Geometric (Unwrap - m)');
grid on; ylabel('Im(s_R)');
legend('Location', 'best');
title('Matching Manual Resets using Accumulated Phase');

subplot(2,1,2)
plot(ImagEigManual - ImagGeometricSawtooth, 'r', 'DisplayName', 'Delta Error');
grid on; ylabel('Error'); xlabel('Index');
legend show;

print('geometric_match_fixed', '-dpng');