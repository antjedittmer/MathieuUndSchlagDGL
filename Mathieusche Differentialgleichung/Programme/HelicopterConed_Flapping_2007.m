% Helicopter Flapping and Dynamic Inflow Analysis
clear; clc; close all;

% Toggle: 1 = 6x6 Cyclic Analysis | 0 = 3x3 Coning Analysis
isSixVec = [0,1];
pos0 = get(0,'defaultFigurePosition');

showDecoupled = 0;

for idxSix = 1:length(isSixVec)
    isSix = isSixVec(idxSix);

    %% 1. Parameters
    Cl_alpha      = 5.7;
    gamma    = 8.159;
    sigma    = 0.07;
    nu0beta     = 1.12;
    CT       = 0.005117;
    lambda_i_bar = sqrt(CT/2);

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
    L_inf = (lambda_i_bar / KI) + C_inf;

    % Following your specific definition for the inflow terms
    C_inf_0 = (Cl_alpha * sigma) / (6 * Km);
    % The term is defined as negative for stability in the A matrix
    L_inf_0 = -(1 / Km) * (4 * lambda_i_bar + (Cl_alpha * sigma / 4));


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

    if isSix == 1
        A_flap_only = A(1:4, 1:4); % Coupled flapping states
        A_inf_only  = A(5:6, 5:6); % Coupled inflow states
        ev_dec = [eig(A_flap_only); eig(A_inf_only)];
    else
        % For 3x3 coning: rows 1-2 are flapping, row 3 is inflow
        ev_dec = [eig(A(1:2, 1:2)); A(3,3)]; % Eventually change this?
    end

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
    disp('Decoupled eigenvalues:'); disp(ev_dec);
    for idx = 1:size(V, 2)
        fprintf('\nMode %d Eigenvector (Mag | Phase Deg):\n', idx);
        disp([abs(V_norm(:,idx)), angle(V_norm(:,idx))*180/pi]);
    end

    % --- Section 3 Expansion: Detailed Eigenvalue Map ---
    figure(idxSix*10);
    %set(gcf,'position',[pos0(1:2) 1.5 *pos0(3), pos0(4)])
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
            k, real(ev(k)), imag(ev(k)), ev_mag, ev_phase); %#ok<*SAGROW>
    end

    grid on; %xline(0, 'k-', 'LineWidth', 2); yline(0, 'k-', 'LineWidth', 2);
    xlabel('Real Part (Damping \sigma)');
    ylabel('Imaginary Part (Frequency \omega)');

    % title('System Eigenvalues: Stability & Frequency Analysis');

    % Capture the handle for the decoupled plot
    h_dec = plot(real(ev_dec), imag(ev_dec), 'kx', 'MarkerSize', 10, 'LineWidth', 1.5);

    % 2. Combine the handles and the labels into one legend call
    % [h; h_dec] merges the 6 mode handles with the 1 decoupled handle
    % [ev_labels, {'Decoupled Eig. Val.'}] merges the two cell arrays
    legend([h; h_dec], [ev_labels, {'Decoupled Eigenvalues'}], ...
        'Location', 'northoutside', 'FontSize', 9, 'FontName', 'Consolas');

    %axis tight;
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
        u = ones(size(t))*deg2rad;
        sys_coned = ss(A, B, eye(3), 0);
        [y, t_out] = lsim(sys_coned, u, t);
    end


    % Create a decoupled version of A by zeroing out coupling blocks
    A_decoupled = A;
    if isSix == 1
        % 6x6 Cyclic Analysis: States [beta_dot_1c/s; beta_1c/s; lambda_1c/s]
        A_decoupled(1:4, 5:6) = 0; % Remove inflow effect on flapping [cite: 57-89]
        A_decoupled(5:6, 1:4) = 0; % Remove flapping effect on inflow [cite: 57-89]

        % Construct the decoupled system for simulation
        sys_dec = ss(A_decoupled, B, eye(size(A, 1)), 0);

        y_dec1 = lsim(sys_dec, u1, t);
        y_dec2 = lsim(sys_dec, u2, t);

    else
        % 3x3 Coning Analysis: States [beta_dot; beta; lambda_i0]

        % Row 1 & 2 are flapping, Row 3 is inflow
        A_decoupled(1:2, 3) = 0;   % Remove inflow effect on flapping
        A_decoupled(3, 1:2) = 0;   % Remove flapping effect on inflow
        % Construct the decoupled system for simulation
        sys_dec = ss(A_decoupled, B, eye(size(A, 1)), 0);
        [y_dec] = lsim(sys_dec, u, t);

        % C_inf_0/(- L_inf_0) * pi/180 - y(end,3) % for debugging
        lambda_inf = C_inf_0/(- L_inf_0) * pi/180;
        u_cor = u - g6* lambda_inf; % Move lambda steady state into input
        [y_dec_cor] = lsim(sys_dec, u_cor, t);
    end



    %% 5. Plotting (States in Degrees)
    lw = 1;
    grey_color = 0.5*ones(1,3);
    if isSix == 1
        % --- 6x6 Cyclic Results: Flapping ---

        figure(2); clf;
        % --- Subplot 1: Theta_C Step ---
        subplot(2,1,1);
        plot(t, y1(:,3)*180/pi, 'b', 'DisplayName', 'Coupled \beta_{1C}','LineWidth', lw); hold on;
        plot(t, y1(:,4)*180/pi, 'r', 'DisplayName', 'Coupled \beta_{1S}','LineWidth', lw);
        if showDecoupled == 1
            plot(t, y_dec1(:,3)*180/pi, 'k:', 'LineWidth', lw+0.5, 'DisplayName', 'Decoupled \beta_{1C}');
            plot(t, y_dec1(:,4)*180/pi, '--', 'Color', grey_color, 'LineWidth', lw+0.5, 'DisplayName', 'Decoupled \beta_{1S}');
        end
        ylabel('Flapping [deg]'); title('Step Response: \theta_C = 1^\circ');
        grid on; legend('Location','east');

        % --- Subplot 2: Theta_S Step ---
        subplot(2,1,2);
        plot(t, y2(:,3)*180/pi, 'b', 'DisplayName', 'Coupled \beta_{1C}','LineWidth', lw); hold on;
        plot(t, y2(:,4)*180/pi, 'r', 'DisplayName', 'Coupled \beta_{1S}', 'LineWidth', lw);
        if showDecoupled == 1
            plot(t, y_dec2(:,3)*180/pi, 'k:', 'LineWidth', lw+0.5, 'DisplayName', 'Decoupled \beta_{1C}');
            plot(t, y_dec2(:,4)*180/pi, '--', 'Color', grey_color, 'LineWidth', lw+0.5, 'DisplayName', 'Decoupled \beta_{1S}');
        end
        ylabel('Flapping [deg]'); xlabel('Time [s]'); title('Step Response: \theta_S = 1^\circ');
        grid on; legend('Location','east');

        figure(3); clf;

        % --- Subplot 1: Theta_C Step ---
        subplot(2,1,1);
        % Note: Removed *180/pi as inflow is typically non-dimensional [-]
        plot(t, y1(:,5), 'b', 'DisplayName', 'Coupled \lambda_{1C}', 'LineWidth', lw); hold on;
        plot(t, y1(:,6), 'r', 'DisplayName', 'Coupled \lambda_{1S}', 'LineWidth', lw);
        if showDecoupled == 1
            plot(t, y_dec1(:,5), 'k:', 'LineWidth', lw+0.5, 'DisplayName', 'Decoupled \lambda_{1C}');
            plot(t, y_dec1(:,6), '--', 'Color', grey_color, 'LineWidth', lw+0.5, 'DisplayName', 'Decoupled \lambda_{1S}');
        end
        ylabel('Inflow [-]'); title('Step Response: \theta_C = 1^\circ');
        grid on; legend('Location','Southeast');

        % --- Subplot 2: Theta_S Step ---
        subplot(2,1,2);
        plot(t, y2(:,5), 'b', 'DisplayName', 'Coupled \lambda_{1C}', 'LineWidth', lw); hold on;
        plot(t, y2(:,6), 'r', 'DisplayName', 'Coupled \lambda_{1S}', 'LineWidth', lw);
        if showDecoupled == 1
            plot(t, y_dec2(:,5), 'k:', 'LineWidth', lw+0.5, 'DisplayName', 'Decoupled \lambda_{1C}');
            plot(t, y_dec2(:,6), '--', 'Color', grey_color, 'LineWidth', lw+0.5, 'DisplayName', 'Decoupled \lambda_{1S}');
        end
        ylabel('Inflow [-]'); xlabel('Time [s]'); title('Step Response: \theta_S = 1^\circ');
        pos0 = get(gca,'YLim');
        set(gca, 'YLim', round(pos0*1000)/1000)
        drawnow;

        grid on; legend('Location','northeast');

    else
        % --- 3x3 Coning Results - SINGLE PLOT MODE ---

        figure(4); clf;

        % --- Subplot 1: Flapping ---
        subplot(2,1,1);
        plot(t_out, y(:,2)*180/pi, 'b', 'LineWidth', lw, ...
            'DisplayName', 'Coupled');
        hold on;
        plot(t_out, y_dec(:,2)*180/pi, 'k:', 'LineWidth', lw+0.5, ...
            'DisplayName', 'Decoupled');
        plot(t_out, y_dec_cor(:,2)*180/pi, '--', 'Color', grey_color, 'LineWidth', lw+0.5, ...
            'DisplayName', ['Decoupled', sprintf(' %2.2f deg', unique(u_cor)*180/pi)]);

        ylabel('\beta_0 [deg]');
        title('Step Response: \theta_0 = 1^\circ');
        grid on;
        legend('Location','best'); % Automatically uses 'DisplayName' labels

        % --- Subplot 2: Inflow ---
        subplot(2,1,2);
        plot(t_out, y(:,3), 'r', 'LineWidth', lw, ...
            'DisplayName', 'Coupled');
        hold on;

        plot(t_out, y_dec(:,3), 'k:', 'LineWidth', lw+0.5, ...
            'DisplayName', 'Decoupled');
        if showDecoupled == 1
            plot(t_out, y_dec_cor(:,3), '--', 'Color', grey_color, 'LineWidth', lw+0.5,...
                'DisplayName', ['Decoupled', sprintf(' %2.2f deg', unique(u_cor)*180/pi)]);
        end

        ylabel('\lambda_{i0} [-]');
        xlabel('Time [s]');
        grid on;
        legend('Location','best');

    end
end