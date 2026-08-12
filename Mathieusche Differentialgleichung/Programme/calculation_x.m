%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mathieu_Floquet_Task.m
%
% Task:
% Compute periodic eigenvector matrix A(t) / V(t) and solution x(t)
% for nu_c = nu_0 = {0.5, 5, 8}, at 63 values t = 0:0.1:2*pi.
% Compare Floquet reconstruction with direct ODE solution.
%
% Based on:
% - Arnold-style transition matrix / monodromy matrix method
% - Peters Floquet interpretation:
%   x(t) = A(t) * exp(Lambda*t) * A0^{-1} * x0
%   with A(t) periodic, A(0)=A0, A(T)=A0
%
% Important:
% This script interprets the task values as nu0 = nuC = 0.5, 5, 8
% in the equation
%   phi_ddot + 2D phi_dot + (nu0^2 + nuC^2*cos(t))*phi = 0
%
% If your professor instead means nu0^2 = nuC^2 = {0.5,5,8},
% then set useSquaredInput = true below.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% User settings
D = 0.15;
nuInputList = [0.5, 5, 8];
useSquaredInput = false;   % false: input means nu0,nuC ; true: input means nu0^2,nuC^2

T = 2*pi;
t0 = 0;
tGrid = 0:0.1:2*pi;        % 63 points
Nz = 2;

% Initial condition for comparison x(t)
% Replace if your task sheet specifies another one
x0 = [1; 0];

% Numerical solver options
opts = odeset('RelTol',1e-10,'AbsTol',1e-12);

% Output containers
AllResults = struct();

%% Loop over the three requested cases
for caseIdx = 1:length(nuInputList)

    nuIn = nuInputList(caseIdx);

    if useSquaredInput
        nu02 = nuIn;
        nuC2 = nuIn;
        nu0  = sqrt(nu02);
        nuC  = sqrt(nuC2);
    else
        nu0  = nuIn;
        nuC  = nuIn;
        nu02 = nu0^2;
        nuC2 = nuC^2;
    end

    fprintf('\n=====================================================\n');
    fprintf('Case %d of %d: nu0 = nuC = %.6g\n', caseIdx, length(nuInputList), nuIn);
    fprintf('Using coefficients nu0^2 = %.6g, nuC^2 = %.6g\n', nu02, nuC2);
    fprintf('=====================================================\n');

    %% Step 1: Compute principal transition matrix Phi(t)
    % Phi(0)=I, built column-by-column by integrating unit initial conditions
    Phi_t = zeros(Nz, Nz, length(tGrid));
    Monodromy = zeros(Nz, Nz);
    I2 = eye(Nz);

    for k = 1:Nz
        solBasis = ode45(@(t,x) MathieuDGL_task(t, x, D, nu02, nuC2), [t0, T], I2(:,k), opts);

        for j = 1:length(tGrid)
            Phi_t(:,k,j) = deval(solBasis, tGrid(j));
        end

        Monodromy(:,k) = deval(solBasis, T);
    end

    %% Step 2: Eigenanalysis of monodromy matrix
    % Monodromy = A0 * exp(Lambda*T) * A0^{-1}
    [A0, Mu] = eig(Monodromy);
    multipliers = diag(Mu);

    % Principal characteristic exponents
    % Peters: lambda_j = (1/T) log(mu_j), branch choice is bookkeeping,
    % but the total product A(t)exp(lambda t) is unique
    lambda = (1/T) * log(multipliers);

    % diagonal exponent matrix
    Lambda = diag(lambda);

    % Alternative matrix-log version
    % R = (1/T)*logm(Monodromy);
    R = (1/T) * logm(Monodromy);

    %% Step 3: Compute periodic eigenvector matrix A(t)
    % Peters formula:
    % A(t) = Phi(t) * A0 * exp(-Lambda*t)
    %
    % Here each column is the periodic Floquet eigenvector for the selected
    % branch of the corresponding exponent.
    A_t = zeros(Nz, Nz, length(tGrid));

    for j = 1:length(tGrid)
        tj = tGrid(j);
        A_t(:,:,j) = Phi_t(:,:,j) * A0 * diag(exp(-lambda*tj));
    end

    %% Step 4: Check periodicity A(T)=A(0)
    A0_from_formula = A_t(:,:,1);
    AT_from_formula = A_t(:,:,end);

    periodicityError = norm(AT_from_formula - A0_from_formula, 'fro');

    %% Step 5: Reconstruct x(t) from Floquet form
    % x(t) = A(t) * exp(Lambda*t) * A0^{-1} * x0
    X_floquet = zeros(Nz, length(tGrid));

    A0inv = inv(A0);
    c0 = A0inv * x0;

    for j = 1:length(tGrid)
        tj = tGrid(j);
        X_floquet(:,j) = A_t(:,:,j) * diag(exp(lambda*tj)) * c0;
    end

    %% Step 6: Direct ODE solution for the same initial condition
    solDirect = ode45(@(t,x) MathieuDGL_task(t, x, D, nu02, nuC2), [t0, T], x0, opts);
    X_direct = deval(solDirect, tGrid);

    %% Step 7: Error between direct and Floquet reconstruction
    X_error = X_floquet - X_direct;
    maxErr = max(vecnorm(X_error, 2, 1));

    %% Step 8: Build modal waveform matrix Q(t)
    % Peters often discusses the modal waveform Q(t)=A(t)*exp(-i*omega*t)
    % using only the imaginary part of lambda.
    %
    % For each mode k:
    % Q_k(t) = A_k(t) * exp(-i*imag(lambda_k)*t)
    %
    % This is useful to see harmonic content without exponential damping/growth.
    Q_t = zeros(Nz, Nz, length(tGrid));

    for j = 1:length(tGrid)
        tj = tGrid(j);
        Q_t(:,:,j) = A_t(:,:,j) * diag(exp(-1i*imag(lambda)*tj));
    end

    %% Step 9: Tables for export / inspection
    % 9a) Time history of x(t)
    TimeTable = table(tGrid(:), ...
                      real(X_floquet(1,:))', imag(X_floquet(1,:))', ...
                      real(X_floquet(2,:))', imag(X_floquet(2,:))', ...
                      X_direct(1,:)', X_direct(2,:)', ...
                      real(X_error(1,:))', imag(X_error(1,:))', ...
                      real(X_error(2,:))', imag(X_error(2,:))', ...
                      'VariableNames', { ...
                      't', ...
                      'xFloquet1_real','xFloquet1_imag', ...
                      'xFloquet2_real','xFloquet2_imag', ...
                      'xDirect1','xDirect2', ...
                      'err1_real','err1_imag','err2_real','err2_imag'});

    % 9b) Flatten A(t) into table
    A11 = squeeze(A_t(1,1,:)); A12 = squeeze(A_t(1,2,:));
    A21 = squeeze(A_t(2,1,:)); A22 = squeeze(A_t(2,2,:));

    ATable = table(tGrid(:), ...
                   real(A11), imag(A11), ...
                   real(A12), imag(A12), ...
                   real(A21), imag(A21), ...
                   real(A22), imag(A22), ...
                   'VariableNames', { ...
                   't', ...
                   'A11_real','A11_imag', ...
                   'A12_real','A12_imag', ...
                   'A21_real','A21_imag', ...
                   'A22_real','A22_imag'});

    % 9c) Modal waveform table Q(t)
    Q11 = squeeze(Q_t(1,1,:)); Q12 = squeeze(Q_t(1,2,:));
    Q21 = squeeze(Q_t(2,1,:)); Q22 = squeeze(Q_t(2,2,:));

    QTable = table(tGrid(:), ...
                   real(Q11), imag(Q11), ...
                   real(Q12), imag(Q12), ...
                   real(Q21), imag(Q21), ...
                   real(Q22), imag(Q22), ...
                   'VariableNames', { ...
                   't', ...
                   'Q11_real','Q11_imag', ...
                   'Q12_real','Q12_imag', ...
                   'Q21_real','Q21_imag', ...
                   'Q22_real','Q22_imag'});

    % 9d) Summary table for eigen information
    EigTable = table((1:2)', ...
                     real(multipliers), imag(multipliers), abs(multipliers), ...
                     real(lambda), imag(lambda), ...
                     'VariableNames', { ...
                     'Mode','Multiplier_real','Multiplier_imag','Multiplier_abs', ...
                     'Exponent_real','Exponent_imag'});

    %% Step 10: Save in result struct
    AllResults(caseIdx).nuInput = nuIn;
    AllResults(caseIdx).nu0 = nu0;
    AllResults(caseIdx).nuC = nuC;
    AllResults(caseIdx).nu02 = nu02;
    AllResults(caseIdx).nuC2 = nuC2;
    AllResults(caseIdx).tGrid = tGrid;
    AllResults(caseIdx).x0 = x0;
    AllResults(caseIdx).Phi_t = Phi_t;
    AllResults(caseIdx).Monodromy = Monodromy;
    AllResults(caseIdx).A0 = A0;
    AllResults(caseIdx).multipliers = multipliers;
    AllResults(caseIdx).lambda = lambda;
    AllResults(caseIdx).Lambda = Lambda;
    AllResults(caseIdx).R = R;
    AllResults(caseIdx).A_t = A_t;
    AllResults(caseIdx).Q_t = Q_t;
    AllResults(caseIdx).X_floquet = X_floquet;
    AllResults(caseIdx).X_direct = X_direct;
    AllResults(caseIdx).X_error = X_error;
    AllResults(caseIdx).periodicityError = periodicityError;
    AllResults(caseIdx).maxErr = maxErr;
    AllResults(caseIdx).TimeTable = TimeTable;
    AllResults(caseIdx).ATable = ATable;
    AllResults(caseIdx).QTable = QTable;
    AllResults(caseIdx).EigTable = EigTable;

    %% Step 11: Print summary
    fprintf('Monodromy matrix M = Phi(T):\n');
    disp(Monodromy);

    fprintf('Characteristic multipliers mu:\n');
    disp(multipliers);

    fprintf('Characteristic exponents lambda = (1/T)log(mu):\n');
    disp(lambda);

    fprintf('||A(T)-A(0)||_F = %.6e\n', periodicityError);
    fprintf('max ||x_Floquet(t)-x_direct(t)||_2 = %.6e\n', maxErr);

    %% Step 12: Plots
    figure('Name', sprintf('Case %d: nu=%.3g', caseIdx, nuIn), 'Color', 'w', ...
           'Position', [100 100 900 900]);

    tiledlayout(4,1);

    % x1 comparison
    nexttile;
    plot(tGrid, real(X_floquet(1,:)), 'b-', 'LineWidth', 1.4); hold on;
    plot(tGrid, X_direct(1,:), 'r--', 'LineWidth', 1.2);
    grid on;
    ylabel('x_1(t)');
    title(sprintf('\\nu_0 = \\nu_C = %.3g, D = %.2f', nuIn, D));
    legend('Floquet reconstruction','Direct ODE','Location','best');

    % x2 comparison
    nexttile;
    plot(tGrid, real(X_floquet(2,:)), 'b-', 'LineWidth', 1.4); hold on;
    plot(tGrid, X_direct(2,:), 'r--', 'LineWidth', 1.2);
    grid on;
    ylabel('x_2(t)');
    legend('Floquet reconstruction','Direct ODE','Location','best');

    % error
    nexttile;
    plot(tGrid, real(X_error(1,:)), 'k-', 'LineWidth', 1.2); hold on;
    plot(tGrid, real(X_error(2,:)), 'm--', 'LineWidth', 1.2);
    grid on;
    ylabel('error');
    legend('e_1','e_2','Location','best');

    % periodic eigenvectors A(t): first-row components as example
    nexttile;
    plot(tGrid, real(A11), 'b-', 'LineWidth', 1.2); hold on;
    plot(tGrid, real(A12), 'r--', 'LineWidth', 1.2);
    grid on;
    xlabel('t');
    ylabel('Re(A_{1j}(t))');
    legend('Re(A_{11})','Re(A_{12})','Location','best');

end

%% Save everything
save('Mathieu_Floquet_Task_Results.mat', 'AllResults');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Local function: first-order Mathieu system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dx = MathieuDGL_task(t, x, D, nu02, nuC2)
    dx = zeros(2,1);
    dx(1) = x(2);
    dx(2) = -(nu02 + nuC2*cos(t))*x(1) - 2*D*x(2);
end