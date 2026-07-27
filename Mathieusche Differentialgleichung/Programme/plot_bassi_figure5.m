clear; clc; close all;

% Approximate curve data from the figure
Omega = linspace(200, 230, 400);
wc = 5.96 + 0.0204*(Omega - 200) + 0.00042*(Omega - 200).^2;

% Approximate scatter data
Omega_black = [214.2 215.0 215.8 216.6 217.4 218.2 219.0 219.8 220.6 221.4];
wc_black    = [6.43  6.45  6.48  6.51  6.54  6.57  6.60  6.63  6.66  6.69];

Omega_green = [214.5 215.4 216.3 217.2 218.1 219.0 219.9 220.8 221.7];
wc_green    = [6.44  6.47  6.50  6.53  6.56  6.59  6.62  6.65  6.68];

Omega_blue  = [214.7 215.6 216.5 217.4 218.3 219.2 220.1 221.0 221.9];
wc_blue     = [6.43  6.46  6.49  6.52  6.55  6.58  6.61  6.64  6.67];

% Time test points
Omega_red   = 215.2; wc_red   = 6.50;
Omega_cyan  = 220; wc_cyan  = 6.69;

figure('Color','w');
hold on; box on; grid on;

% Plot and save handles
h1 = plot(Omega, wc, 'Color', [0.65 0.1 0.2], 'LineWidth', 2.5);
h2 = plot(Omega_black, wc_black, 'k+', 'MarkerSize', 10, 'LineWidth', 2);
h3 = plot(Omega_green, wc_green, 'go', 'MarkerSize', 9, 'LineWidth', 2);
h4 = plot(Omega_blue,  wc_blue,  'bs', 'MarkerSize', 9, 'LineWidth', 2);
h5 = xline(210, '--', 'Color', [0.95 0.35 0.35], 'LineWidth', 2.5);
h6 = plot(Omega_red,  wc_red,  'r*', 'MarkerSize', 12, 'LineWidth', 2);
h7 = plot(Omega_cyan, wc_cyan, 'g*', 'MarkerSize', 12, 'LineWidth', 2);

% Axes formatting
xlim([200 230]);
ylim([5.5 7.5]);
xlabel('\Omega [rpm]', 'FontSize', 16, 'Interpreter', 'tex');
ylabel('\omega, \omega_c [rad/s]', 'FontSize', 16, 'Interpreter', 'tex');

set(gca, 'FontSize', 14, 'LineWidth', 1.0, 'TickDir', 'in');
xticks(200:10:230);
yticks(5.5:0.5:7.5);

% Correct legend: use handles in desired order
legend([h1 h2 h3 h4 h5 h6 h7], ...
       {'Unstable mode 5', ...
        '\omega_c, \delta = 0.05', ...
        '\omega_c, \delta = 0.1', ...
        '\omega_c, \delta = 0.15', ...
        'Left instability bound', ...
        'Time test \omega_c', ...
        'Time test \omega_c'}, ...
       'Location', 'eastoutside', ...
       'FontSize', 12);

% Plot area similar to screenshot
set(gca, 'Position', [0.12 0.12 0.56 0.74]);