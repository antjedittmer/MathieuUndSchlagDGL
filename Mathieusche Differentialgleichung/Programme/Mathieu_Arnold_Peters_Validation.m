%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mathieu_Floquet_Verification.m
%
% Task: Reconstruct x(t) (2x63) via Peters V(t) and Arnold s_R, 
% then validate against direct ODE integration.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;

%% 1. User Settings & Inputs
D = 0.15;                       % Damping ratio
nuInputList = [0.5, 5, 8];      % Frequency parameters (nu_0 = nu_C)
mFactorList = [0.5, 2, 2.5];    % Arnold addition factors for the shift
T = 2*pi;                       % Parametric period (T = 2*pi)
t0 = 0;
tGrid = 0:0.1:2*pi;             % 63 distinct timesteps (1x63 vector)
Nz = 2;                         % Number of state dimensions

% Complete list of required initial conditions
x0_cases = {[1; 0], [0; 1]};

opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
AllResults = struct();
structIdx = 1;

%% 2. Execution Loop across Parameter Cases
for caseIdx = 1:length(nuInputList)
    nuIn = nuInputList(caseIdx);
    m_factor = mFactorList(caseIdx);
    
    nu02 = nuIn^2;
    nuC2 = nuIn^2;
    
    fprintf('\n=====================================================\n');
    fprintf('Case %d: nu_0 = nu_C = %.1f | Arnold m = %.1f\n', caseIdx, nuIn, m_factor);
    fprintf('=====================================================\n');
    
    %% Step A: Compute the fundamental transition matrix Phi(t)
    Phi_t = zeros(Nz, Nz, length(tGrid));
    Monodromy = zeros(Nz, Nz);
    I2 = eye(Nz);
    
    for k = 1:Nz
        % Integrate standard basis unit vectors over the full period
        solBasis = ode45(@(t,x) MathieuDGL_task(t, x, D, nu02, nuC2), [t0, T], I2(:,k), opts);
        for j = 1:length(tGrid)
            Phi_t(:,k,j) = deval(solBasis, tGrid(j));
        end
        Monodromy(:,k) = deval(solBasis, T); 
    end
    
    %% Step B: Eigenanalysis & Arnold Exponent Matrix Calculation
    [V0, Mu] = eig(Monodromy);
    multipliers = diag(Mu);
    
    % Raw characteristic exponents via the principal branch logarithm
    s_raw = (1/T) * log(multipliers);
    
    % Formulate the characteristic exponents according to Arnold
    % Shifts are applied to conjugate pairs to trace continuous paths
    s_R = zeros(Nz, 1);
    s_R(1) = real(s_raw(1)) + 1i * (imag(s_raw(1)) - m_factor);
    s_R(2) = real(s_raw(2)) + 1i * (imag(s_raw(2)) + m_factor);
    
    fprintf('  Arnold Exponents s_R:\n');
    fprintf('    s_R1 = %.4f + %.4fi\n', real(s_R(1)), imag(s_R(1)));
    fprintf('    s_R2 = %.4f + %.4fi\n', real(s_R(2)), imag(s_R(2)));
    
    %% Step C: Compute Time-Varying Eigenvector Matrix V(t) (Peters' A(t))
    % To guarantee invariance under branch choices, V(t) shifts in opposition:
    % V(t) = Phi(t) * V(0) * diag(exp(-s_R * t))
    V_t = zeros(Nz, Nz, length(tGrid));
    for j = 1:length(tGrid)
        tj = tGrid(j);
        V_t(:,:,j) = Phi_t(:,:,j) * V0 * diag([exp(-s_R(1)*tj), exp(-s_R(2)*tj)]);
    end
    
    % Verify structural periodicity of Peters' matrix
    periodicityError = norm(V_t(:,:,end) - V_t(:,:,1), 'fro');
    fprintf('  ||V(T) - V(0)||_F Periodicity Error = %.6e\n', periodicityError);
    
    %% Step D: Trajectory Verification for Initial Conditions
    for initIdx = 1:length(x0_cases)
        x0 = x0_cases{initIdx};
        
        % Preallocate target variables (2 x 63 elements)
        X_floquet = zeros(Nz, length(tGrid));
        V0inv = inv(V0);
        
        for j = 1:length(tGrid)
            tj = tGrid(j);
            % Structured Reconstruction Equation
            exp_matrix = diag([exp(s_R(1)*tj), exp(s_R(2)*tj)]);
            X_floquet(:,j) = V_t(:,:,j) * exp_matrix * V0inv * x0;
        end
        % Filter tiny machine-epsilon imaginary noise from real physical space
        X_floquet = real(X_floquet); 
        
        % Direct ground-truth solution via ODE solver
        solDirect = ode45(@(t,x) MathieuDGL_task(t, x, D, nu02, nuC2), [t0, T], x0, opts);
        X_direct = deval(solDirect, tGrid);
        
        % Calculate global trajectory error
        X_error = X_floquet - X_direct;
        maxErr = max(vecnorm(X_error, 2, 1));
        fprintf('  Initial Condition [%d; %d] -> Max Abs Error = %.6e\n', x0(1), x0(2), maxErr);
        
        %% Step E: Save to global struct matrix
        AllResults(structIdx).nu = nuIn;
        AllResults(structIdx).m = m_factor;
        AllResults(structIdx).x0 = x0;
        AllResults(structIdx).s_R = s_R;
        AllResults(structIdx).V_t = V_t;
        AllResults(structIdx).X_floquet = X_floquet; % Matrix of 2x63 elements
        AllResults(structIdx).X_direct = X_direct;   % Matrix of 2x63 elements
        AllResults(structIdx).maxErr = maxErr;
        
        structIdx = structIdx + 1;
    end
end

% Save verified results matrix
save('Mathieu_Floquet_Verification_Results.mat', 'AllResults');
disp('=====================================================\n');
disp('Verification completed successfully! Data stored in workspace.');

%% Local Function: First-Order Mathieu State Space System
function dx = MathieuDGL_task(t, x, D, nu02, nuC2)
    dx = zeros(2,1);
    dx(1) = x(2); % Position derivative = velocity
    dx(2) = -(nu02 + nuC2*cos(t))*x(1) - 2*D*x(2); % Parametric Mathieu oscillator equation
end