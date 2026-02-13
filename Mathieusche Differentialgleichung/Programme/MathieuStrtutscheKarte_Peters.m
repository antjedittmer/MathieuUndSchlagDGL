%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stability Analysis (Strutt Chart) using Floquet Method.
%
% Arnold-style presentation with Peters-style frequency tracking.
% phi'' + 2D*phi' + (nu_0^2 + nu_C^2*cos(psi))*phi = 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;

% --- Setup and Parameters ---

loadMat = 1; % Load mat-file if results with the same D are already available
SW = 0.05;    % Step width for nu_0^2 and nu_C^2 sweep
unt0 = 0;    % Lower bound for nu_0^2 (x-axis)
untC = 0;    % Lower bound for nu_C^2 (y-axis)
ob0  = 9;    % Upper bound for nu_0^2
obC  = 9;    % Upper bound for nu_C^2

% --- Setup for Figure Saving ---
fDir = 'figureFolder'; % Folder for figures
if ~isdir(fDir) %#ok<ISDIR>
    mkdir(fDir)
end
fDirPeters = fullfile(fDir,'figureFolder_Arnold_Peters'); % Subfolder specific to Peters' plots
if ~isdir(fDirPeters) %#ok<ISDIR>
    mkdir(fDirPeters)
end
dDir = 'dataFolder'; % Folder for figures
if ~isdir(dDir) %#ok<ISDIR>
    mkdir(dDir)
end
dDir1 = fullfile(dDir,'dataFolder_Arnold_Peters');
if ~isfolder(dDir1)
    mkdir(dDir1)
end
excelDir =fullfile(dDir,'dataFolder_Arnold_Peters_Excel');
if ~isfolder(excelDir)
    mkdir(excelDir);
end

Nz   = 2;      % Number of equations in the DGL system (2nd-order ODE)
DVec = 0.15;   % Damping coefficient D
t0    = 0.0;
Omega = 1;          % Normalized excitation frequency
T     = 2*pi/Omega; % Period of excitation

for dIdx = 1:length(DVec)
    D = DVec(dIdx);

    % Initial conditions (identity)
    Diagonal = diag(ones(Nz,1));

    % MAT filename
    matName  = [strrep(sprintf('STRUTTscheKarte_Arnold_Peters_D%2.1e_SW%2.1e',D,SW),'.','dot'),'.mat'];
    fileName = fullfile(dDir,matName);

    % Grid
    nu02_vals = unt0:SW:ob0;
    nuC2_vals = untC:SW:obC;
    lenNu02   = length(nu02_vals);
    lenNuC2   = length(nuC2_vals);
    lenNu     = lenNu02 * lenNuC2;
    lenNuDiag = min(lenNu02,lenNuC2);

    buffer.Pos = 0;
    buffer.Neg = 0;

    % Storage: nu02, nuC2, Re_s1, Re_s2, Im_s1/Omega, Im_s2/Omega
    CharEx         = zeros(lenNuDiag, 6);
    plotwertstabil = zeros(lenNu,3);
    lidx = 1;
    oidx = 1;
    m_direct = zeros(lenNuDiag + 1,1);
    m_new = zeros(lenNuDiag + 1,1);
    cnt = 1;

    if exist(fileName,'file') == 2 && loadMat == 1
        load(fileName,'CharEx','plotwertstabil');
    end
    if ~exist('CharEx','var') || size(CharEx,2) ~= 8
        disp(['Calculating Strutt Chart (Arnold/Peters Style) for D = ', num2str(D), '...']);
        for nu_02 = nu02_vals
            omega0 = sqrt(nu_02); % undisturbed natural frequency
            for nu_C2 = nuC2_vals
                options = odeset('RelTol',1e-10,'AbsTol',1e-12);

                % Monodromy matrix
                Monodromie = zeros(Nz);
                for k = 1:Nz
                    sol = ode45(@(psi,x)MathieuDGLsubfun(psi,x,D,nu_02,nu_C2), [t0,T], Diagonal(:,k), options);
                    Monodromie(:,k)        = deval(sol,T);


                    % 3. Extract the winding number from the FIRST basis vector (k=1)
                    if k == 1 && abs(nu_C2 - nu_02) < SW/2
                        % Sample the trajectory over the period
                        psi_eval = linspace(t0, T, 500);
                        z_path = deval(sol, psi_eval); % z_path(1,:) is x_dot, z_path(2,:) is x

                        % Calculate polar angle in the (x, x_dot) plane
                        % Note: atan2(y, x) -> atan2(x_dot, x)
                        % This corresponds to atan2(z_path(1,:), z_path(2,:))
                        theta = unwrap(atan2(z_path(1,:), z_path(2,:)));

                        % Total accumulated phase change
                        total_delta_theta = abs(theta(end) - theta(1));

                        % m is the count of half-rotations (0.5 per 180 degrees)
                        m_direct(oidx) = round(total_delta_theta / pi) * 0.5;
                    end

                end

                % Characteristic multipliers
                eP = eig(Monodromie);

                % --- Characteristic exponents (diagonal only) ---
                if abs(nu_C2 - nu_02) < SW/2
                    Eig.Real = 1/T * log(abs(eP));
                    Eig.Imag = 1/T * atan(imag(eP)./real(eP));

                    % continuity correction
                    [Eig,buffer] = correctImagValues(Eig,buffer);

                    % Peters-style physical frequencies (here just stacked)
                    PhysFreq = [Eig.ImagCorrected; Eig.ImagCorrectedNeg];

                    % sort by real part
                    Eig_Re = Eig.Real;
                    [~, idx_sort] = sort(Eig_Re);

                    % store: nu02, nuC2, Re_s1, Re_s2, Im_s1, Im_s2
                    CharEx(oidx,:) = [nu_02, nu_C2, Eig_Re(idx_sort)', (PhysFreq)'];
                    oidx = oidx + 1;

                    % Check if we are in a "bubble": Real parts are not identical
                    isBubble = abs(Eig_Re(1) - Eig_Re(2)) > eps;

                    if isBubble
                        % If "bubble", allow m to follow the ODE winding number
                        m_new(oidx) = m_direct(oidx-1);
                        cnt = cnt +1;
                    else
                        % If not inside a bubble, lock m to its entry value
                        % This prevents the 'early jump' during the resonance
                        m_new(oidx) = m_new(oidx - 1);
                    end
                end

                % --- Stability map ---
                if max(abs(eP)) < 1
                    plotwertstabil(lidx,:) = [nu_02,nu_C2,1];
                end
                lidx = lidx + 1;
            end
        end
        CharEx = [CharEx, m_direct(2:end), m_new(2:end)]; %#ok<AGROW>
        % remove zero rows
        plotwertstabil = plotwertstabil(any(plotwertstabil,2),:);
        save(fileName,'CharEx','plotwertstabil');
    end

    % --- Excel export ---
    try
        CharExTable = array2table(CharEx,'VariableNames',...
            {'nu02','nuC2','Re_s1','Re_s2', 'Freq_s1_norm', 'Freq_s2_norm','m_direct','m_new'});
        excelfilename  = strrep(fileName,'.mat','_CharExPhys.xlsx');
        excelfilename1 = strrep(excelfilename,dDir,excelDir);
        writetable(CharExTable,excelfilename1)
    catch ME
        warning(['Error during Excel export: ', ME.message]);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Graphical Representation (Arnold/Peters + test function)
    cl = lines;
    fs = 10.5;
    hf = figure(dIdx);
    hf.Position = [10 10 600 700];

    % Diagonal parameter
    xachse = CharEx(:,1);

    % -------------------- 1) Strutt Chart ---------------------------
    h(1) = subplot(3,1,1);
    scatter(plotwertstabil(:,1),plotwertstabil(:,2),5,'filled');
    title({'Strutt Chart for $\ddot{\phi} + 2D \dot{\phi} + (\nu^2_0 + \nu^2_C \cos(\psi))\phi = 0$'},...
        'interpreter','latex','FontSize', fs+2);
    ylabel('Parameter $\nu_C^2$ $\rm{[-]}$','interpreter','latex','FontSize', fs+2);
    grid on;

    % -------------------- 2) Real parts -----------------------------
    h(2) = subplot(3,1,2);
    plot(xachse,CharEx(:,3),'LineWidth', 1.5, 'Color', cl(1,:), 'DisplayName', '$\sigma_1$');
    hold on;
    plot(xachse,CharEx(:,4),'LineWidth', 1.5, 'Color', cl(2,:), 'DisplayName', '$\sigma_2$');
    yline(0, 'k--');
    title('Real Part (Damping) $\sigma = Re(s) = \frac{1}{T}\ln(|\mu|)$','interpreter','latex','FontSize', fs+2);
    ylabel('Damping Exponent $\sigma$ $\rm{[-]}$','interpreter','latex','FontSize', fs+2);
    grid on;
    legend('Location','NorthEast', 'Interpreter','latex');

    % -------------------- 3) Imaginary part: test function ---------
    h(3) = subplot(3,1,3);

    % Build candidate lines and omega from CharExTable
    Im_AllLine = CharExTable.Freq_s1_norm + (0:0.5:2.5).*ones(size(CharExTable.Freq_s1_norm));
    omega      = sqrt(CharExTable.nu02);

    % closest-from-below + monotone selection
    Im_m   = NaN(size(omega));
    Im_m_d = NaN(size(omega));
    idxmin = NaN(size(omega));

    for idx = 1:length(omega)
        diff = Im_AllLine(idx,:) - omega(idx);   % no abs
        mask = diff <= 0.1;                        % only from below

        [Im_m(idx), idx_loc] = max(Im_AllLine(idx,mask));
        idx_all              = find(mask);
        idxmin(idx)          = idx_all(idx_loc);
        Im_m_d(idx)          = omega(idx) - Im_m(idx);
    end

    plot(xachse, Im_m, 'LineWidth', 1.0, ...
        'Color', cl(1,:), 'DisplayName', '$\omega$ nearest smaller $\sqrt{\nu_c^2}$');

    hold on; grid on;

    plot(xachse,CharExTable.m_new + CharExTable.Freq_s1_norm,'LineWidth', 1.3, 'Color', cl(2,:),...
        'LineStyle',':','DisplayName','$\omega$ with $m$ from basis vector')
    plot(xachse, omega, 'k--', 'LineWidth', 1.0, 'DisplayName', '$\sqrt{\nu_c^2}$');
    title('Imaginary Part (test function)','interpreter','latex','FontSize', fs+2)
    xlabel('Parameter $\nu_0^2$ $\rm{[-]}$','interpreter','latex','FontSize', fs+2);
    ylabel('Im$(s_R)$ $\rm{[-]}$','interpreter','latex','FontSize', fs+2);

    legend('Location','SouthEast', 'Interpreter','latex');

    linkaxes(h,'x')
    for idxH = 1:length(h)
        set(h(idxH),'TickLabelInterpreter','Latex','FontSize',fs)
    end

    pngname = fullfile(fDirPeters, strrep(matName, '.mat', '_Exponents'));
    print(pngname, '-dpng')
end

% =========================================================================
% AUXILIARY FUNCTIONS
% =========================================================================

function PhysFreq = petersPhysicalFrequency(Eig_Im_Raw, ~, omega0, Omega)
% petersPhysicalFrequency: picks integer k so Im(s)+k*Omega is closest to omega0.
k_range  = -5:5;
PhysFreq = zeros(size(Eig_Im_Raw));
for i = 1:length(Eig_Im_Raw)
    nu_imag_raw      = Eig_Im_Raw(i);
    test_frequencies = nu_imag_raw + k_range * Omega;
    target_freq      = omega0;
    [~, k_best_idx]  = min(abs(test_frequencies - target_freq));
    PhysFreq(i)      = test_frequencies(k_best_idx);
end
end

function dxdpsi = MathieuDGLsubfun(psi, x, D, nu_02, nu_C2)
% phi'' + 2D*phi' + (nu_0^2 + nu_C^2*cos(psi))*phi = 0
phi     = x(1);
phi_dot = x(2);
K_psi   = nu_02 + nu_C2 * cos(psi);
phi_ddot = -2 * D * phi_dot - K_psi * phi;
dxdpsi  = [phi_dot; phi_ddot];
end

function  [Eig, buffer] = correctImagValues(Eig, buffer)
% correctImagValues: makes Imag parts continuous in parameter sweep.
if nargin~= 2
    error('Two inputs expected: Eig struct and buffer struct.');
end
if max(size(Eig.Imag))~= 2 || min(size(Eig.Imag))~= 1
    error('Eig.Imag must contain the pair of imaginary parts.');
end

Eig.ImagSort = sort(Eig.Imag);
tmp    = Eig.ImagSort(2);
tmpNeg = Eig.ImagSort(1);

if ~isfield(buffer, 'Pos'), buffer.Pos = 0; end
if ~isfield(buffer, 'Neg'), buffer.Neg = 0; end

if buffer.Pos <= tmp || (abs(tmp) < 10^-5)
    Eig.ImagCorrected    = tmp;
    Eig.ImagCorrectedNeg = tmpNeg;
    buffer.Pos = 0;
    buffer.Neg = 0;
else
    Eig.ImagCorrected    = 2*buffer.Pos + tmpNeg;
    Eig.ImagCorrectedNeg = 2*buffer.Neg + tmp;
end

buffer.Pos = max(tmp,buffer.Pos);
buffer.Neg = min(tmpNeg,buffer.Neg);
end
