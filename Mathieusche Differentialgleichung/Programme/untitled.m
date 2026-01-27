%% Comparison: Arnold vs Peters imaginary parts (as per supervisor)
clc; close all;

% === Folders (match your scripts) ===
dDirA = fullfile('dataFolder', 'dataFolder_Arnold_Classic_Symmetric_test');
dDirP = fullfile('dataFolder', 'dataFolder_Arnold_Peters');  % Adjust if different

% === Load Arnold (trusted arctan Im(s)) ===
matA = 'STRUTTscheKarte_D1dot5e-01_SW1dot0e-01_unt0.mat';
S_A  = load(fullfile(dDirA, matA));
nu_A = S_A.CharEx(:,1);                    % nu_02 (diagonal)
ImA1 = S_A.CharEx(:,7);                    % ImagEig1 (corrected + n)
ImA2 = S_A.CharEx(:,8);                    % ImagEig2 (corrected + n)

% === Load Peters (test function Im_m) ===
% NOTE: Peters script doesn't save Im_m, so you need to either:
% Option 1: Run Peters script first and capture Im_m from workspace
% Option 2: Rebuild it here (code below)
matP = 'STRUTTscheKarte_Arnold_Peters_D1dot5e-01_SW1dot0e-01.mat';
S_P  = load(fullfile(dDirP, matP));
nu_P = S_P.CharEx(:,1);
Freq_s1 = S_P.CharEx(:,5);                 % Freq_s1_norm (Peters physical freq)

% Rebuild Im_m exactly as in Peters bottom subplot:
m_vec = 0:0.5:2.5;
Im_AllLine = Freq_s1 + m_vec .* ones(size(Freq_s1));
omega0     = sqrt(nu_P);
Im_m = NaN(size(omega0));
Im_m_d = NaN(size(omega0));
idxmin = NaN(size(omega0));
for idx = 1:length(omega0)
    diff = Im_AllLine(idx,:) - omega0(idx);
    mask = diff <= 0;  % closest from below
    if any(mask)
        [Im_m(idx), idx_loc] = max(Im_AllLine(idx,mask));
        idx_all = find(mask);
        idxmin(idx) = idx_all(idx_loc);
        Im_m_d(idx) = omega0(idx) - Im_m(idx);
    else
        [Im_m_d(idx), idxmin(idx)] = min(Im_AllLine(idx,:) - omega0(idx));
        Im_m(idx) = Im_AllLine(idx,idxmin(idx));
    end
end

% === Plot comparison ===
figure('Position',[100 100 800 500]); hold on;
plot(nu_A, ImA1, 'b-', 'LineWidth',2, 'DisplayName','Arnold $\operatorname{Im}(s_{R1})$');
plot(nu_A, ImA2, 'c-', 'LineWidth',2, 'DisplayName','Arnold $\operatorname{Im}(s_{R2})$');
plot(nu_P, Im_m, 'r--', 'LineWidth',2, 'DisplayName','Peters test $\operatorname{Im}_m$');
plot(nu_P, omega0, 'k:', 'LineWidth',1.5, 'DisplayName','$\sqrt{\nu_0^2}$');
xlabel('$\nu_0^2$ (diagonal $\nu_C^2 = \nu_0^2$)','Interpreter','latex','FontSize',14);
ylabel('Imaginary part $\operatorname{Im}(s_R)$','Interpreter','latex','FontSize',14);
title('Arnold (trusted) vs Peters test function $\operatorname{Im}(s)$','Interpreter','latex','FontSize',16);
legend('Location','best','Interpreter','latex');
grid on; ylim([0 3]);

% === Check difference (find M?) ===
figure; 
plot(nu_P, Im_m - ImA1(1:length(nu_P)), 'g-', 'LineWidth',2);
hold on; plot(nu_P, Im_m - ImA2(1:length(nu_P)), 'm-', 'LineWidth',2);
yline(0,'k--'); grid on;
xlabel('$\nu_0^2$'); ylabel('Peters Im_m - Arnold Im(s_R)');
title('Difference: constant M?'); legend('vs ImA1','vs ImA2');

fprintf('Check plots: look for constant offset in bubbles/outside bubbles\n');
fprintf('If Im_m ≈ ImA + M, then M is your correction factor.\n');
