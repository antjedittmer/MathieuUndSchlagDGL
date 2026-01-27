
clc; clear; close all;
excelfilename ='STRUTTscheKarte_D1dot5e-01_SW5dot0e-02_unt0CharExAll.xlsx';
foldername = fullfile(fileparts(pwd),'dataFolder','dataFolder_Arnold_Excel_Classic_Symmetric_test');
% load('CharExTableRaw.mat','CharExTable');
CharExTable = readtable(fullfile(foldername,excelfilename));
Im_AllLine = CharExTable.ImagEigNoadditionfactor  + (0:0.5:2.5).*ones(size(CharExTable.ImagEigNoadditionfactor ));
omega      = sqrt(CharExTable.nu02);

figure; hold on;
plot(Im_AllLine);
plot(omega,'k','LineWidth',1.5);
hold off;

Im_m   = NaN(size(omega));
Im_m_d = NaN(size(omega));
idxmin = NaN(size(omega));

for idx = 1:length(omega)
    diff = Im_AllLine(idx,:) - omega(idx);     % no abs
    mask = diff <= 0.1;                          % only from below

    [Im_m(idx), idx_loc] = max(Im_AllLine(idx,mask));
    idx_all              = find(mask);
    idxmin(idx)          = idx_all(idx_loc);
    Im_m_d(idx)          = omega(idx) - Im_m(idx);
end
hold on;
plot(Im_m,'b','LineWidth',1.5);      % blue curve, non-decreasing
plot(omega,'k--','LineWidth',1.5);   % dashed omega
plot(CharExTable.ImagEig2, 'r-.');
xlabel('Index');
ylabel('Value');
grid on;
