clc; clear; close all;

%% Load and assign parameters
load('participation_data_and_M','participation_data', 'm_bubble', ...
    'nu_vals', 'composite_freq')



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
legend('Location', 'best');

% Subplot 4: Winding Number Comparison
subplot(4,1,4)
plot(nu_vals, m_bubble, 'Color', cl(1,:), 'LineWidth', 1.3, 'DisplayName', 'm bubble');
hold on;
plot(nu_vals, m_modpart, '--', 'Color', cl(2,:), 'LineWidth', 1.3, 'DisplayName', 'm Peters');
ylabel('Winding No. m'); grid on; axis tight;
xlabel('Amplification factor \nu');
legend('Location', 'best');