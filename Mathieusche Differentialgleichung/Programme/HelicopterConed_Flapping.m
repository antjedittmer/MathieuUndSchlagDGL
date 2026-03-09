% Helicopter Flapping and Dynamic Inflow Analysis
clear; clc; close all;

% Toggle: 1 = 6x6 Cyclic Analysis | 0 = 3x3 Coning Analysis
isSixVec = [1,0];
for idxSix = 1:2
    isSix = isSixVec(idxSix);

    %% 1. Parameters
    Cl_alpha      = 5.7;
    gamma    = 8.159;
    sigma    = 0.07;
    nu0beta     = 1.12;
    CT       = 0.005117;
    lambda_i = sqrt(CT/2);

    % Using Option 1 for Inflow Constants
    KI = 16 / (pi * 45);
    Km = 8 / (pi * 3); % Note: Km is not used in this specific A-matrix version

    %% 2. Construct Matrix A and B from Images
    % State vector: [beta_dot_1c; beta_dot_1s; beta_1c; beta_1s; lambda_1s; lambda_1c]
    % Inputs: [theta_c; theta_s]

    % Helper variables for clarity
    g8    = gamma / 8;
    g6 = gamma / 6;
    nuSq1 = nu0beta^2 - 1;
    C_inf = (Cl_alpha * sigma) / (16 * KI);
    L_inf = (lambda_i / KI) + C_inf;

    % Following your specific definition for the inflow terms
    C_inf_0 = (Cl_alpha * sigma) / (6 * Km);
    % The term is defined as negative for stability in the A matrix
    L_inf_0 = -(1 / Km) * (4 * lambda_i + (Cl_alpha * sigma / 4));

    if isSix == 1
        A = [ -g8,    -2,    -nuSq1,  -g8,      0,     -g8;
            2,     -g8,     g8,    -nuSq1,  -g8,      0;
            1,      0,      0,      0,      0,      0;
            0,      1,      0,      0,      0,      0;
            0,     -C_inf,  C_inf,  0,     -L_inf,  0;
            -C_inf,  0,      0,     -C_inf,   0,     -L_inf ];
        B = [ g8,      0;
            0,       g8;
            0,       0;
            0,       0;
            0,       C_inf;
            C_inf,   0 ];
    else
        % Matrix A from your derivation (3x3 coning)
        A = [ -g8,      -nu0beta^2,  -g6;
            1,        0,        0;
            -C_inf_0,  0,       L_inf_0 ];
        % Matrix B from your derivation
        B = [ g8;
            0;
            C_inf_0 ];
    end

    %% 3. Eigenvalues and Normalized Eigenvectors
    [V, D] = eig(A);
    ev = diag(D);

    % Custom Normalization:
    % Max element magnitude = 1, Phase of that element = 0 deg
    V_norm = zeros(size(V));
    for k = 1:size(V, 2)
        [max_val, idx] = max(abs(V(:,k)));
        ref_phase = angle(V(idx, k));
        V_norm(:, k) = (V(:, k) / max_val) * exp(-1i * ref_phase);
    end

    % Display Results to Command Window
    disp('Eigenvalues:'); disp(ev);
    for idx = 1:size(V, 2)
        fprintf('\nMode %d Eigenvector (Mag | Phase Deg):\n', idx);
        disp([abs(V_norm(:,idx)), angle(V_norm(:,idx))*180/pi]);
    end

    % --- Section 3 Expansion: Detailed Eigenvalue Map ---
    figure(1);
    clf; hold on;
    colors_ev = lines; % 6 distinct colors for the 6 modes
    h = zeros(size(V, 2),1);
    for k = 1:size(V, 2)
        % Calculate Magnitude and Phase of the Eigenvalue itself
        ev_mag = abs(ev(k));
        ev_phase = angle(ev(k)) * 180/pi; % Phase in degrees

        % Plot each eigenvalue
        h(k) = plot(real(ev(k)), imag(ev(k)), 'o', 'MarkerSize', 8, ...
            'MarkerFaceColor', colors_ev(k,:), 'Color', 'k', 'LineWidth', 1.5);

        % Construct legend string with Mag and Phase
        ev_labels{k} = sprintf('Mode %d: %.2f + %.2fi (Mag: %.2f, Ph: %.1f°)', ...
            k, real(ev(k)), imag(ev(k)), ev_mag, ev_phase);
    end
    grid on; xline(0, 'k-', 'LineWidth', 2); yline(0, 'k-', 'LineWidth', 2);
    xlabel('Real Part (Damping \sigma)');
    ylabel('Imaginary Part (Frequency \omega)');
    title('System Eigenvalues: Stability & Frequency Analysis');
    legend(h, ev_labels, 'Location', 'bestoutside', 'FontSize', 9, 'FontName', 'Consolas');
    axis padded;
    hold off;

    %% 4. Time Simulation (Step Inputs)
    deg2rad = pi/180;
    t = 0:0.01:25; % Non-dimensional time

    if isSix == 1

        % Inputs must be (length(t) x num_inputs)
        % Case 1: QC = 1 deg, QS = 0
        u1 = [ones(length(t), 1)*deg2rad, zeros(length(t), 1)];
        % Case 2: QC = 0, QS = 1 deg
        u2 = [zeros(length(t), 1), ones(length(t), 1)*deg2rad];

        sys = ss(A, B, eye(size(A, 1)), 0);
        [y1, t1] = lsim(sys, u1, t);
        [y2, t2] = lsim(sys, u2, t);
    else
        % 3x3 coning
        u = ones(size(t));

        sys_coned = ss(A, B, eye(3), 0);
        [y, t_out] = lsim(sys_coned, u, t);
    end


    %% 5. Plotting (States in Degrees)
    if isSix == 1
        % --- 6x6 Cyclic Results---
        figure(2);
        subplot(2,1,1);
        plot(t, y1(:,3)*180/pi, 'b', t, y1(:,4)*180/pi, 'r--');
        ylabel('Flapping [deg]'); title('Step Response: \theta_C = 1^\circ');
        legend('\beta_{1C}', '\beta_{1S}'); grid on;
        subplot(2,1,2);
        plot(t, y2(:,3)*180/pi, 'b', t, y2(:,4)*180/pi, 'r--');
        ylabel('Flapping [deg]'); xlabel('Time [s]');
        title('Step Response: \theta_S = 1^\circ');
        legend('\beta_{1C}', '\beta_{1S}'); grid on;

        % 6x6 Inflow plotting
        figure(3);
        subplot(2,1,1);
        plot(t, y1(:,5)*180/pi, 'b', t, y1(:,6)*180/pi, 'r--');
        ylabel('Inflow [-]'); title('Step Response: \theta_C = 1^\circ');
        legend('\lambda_{1C}', '\lambda_{1S}'); grid on;
        subplot(2,1,2);
        plot(t, y2(:,5)*180/pi, 'b', t, y2(:,6)*180/pi, 'r--');
        ylabel('Inflow [-]'); xlabel('Time [s]');
        title('Step Response: \theta_S = 1^\circ');
        legend('\lambda_{1C}', '\lambda_{1S}'); grid on;

    else
        % --- 3x3 Coning Results - SINGLE PLOT MODE ---
        figure(4); clf;

        subplot(2,1,1);
        plot(t_out, y(:,2)*180/pi, 'b', 'LineWidth', 1.5);  % beta_0 flapping
        ylabel('\beta_0 [deg]');
        title('Step Response: \theta_0 = 1^\circ');
        grid on;

        subplot(2,1,2);
        plot(t_out, y(:,3), 'r', 'LineWidth', 1.5);         % lambda_i0 inflow
        ylabel('\lambda_{i0} [-]'); xlabel('Time [s]');
        grid on;

        % NO legend needed - only 1 signal per subplot (exactly as requested)
    end
end