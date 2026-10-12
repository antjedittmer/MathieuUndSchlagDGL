clc; clear; close all;

%% Lade die Daten und erstelle den Table
matFolder = fullfile(fileparts(pwd),'dataFolder'); % mat Datei im dataFolder
% matFolder = fullfile(pwd);% mat Datei im gleichen Ordner
matName = fullfile(matFolder, 'STRUTTscheKarte_D1dot5e-01_SW1dot0e-01_unt0.mat '); 

load(matName,'Char*','plotwert*', 'nu*');

tableCharEx = array2table(CharEx,'VariableNames',...
    {'nu_02','nu_C2','RealEig1','RealEig2','ImagEig1_0','ImagEig2_0',...
    'ImagEig1','ImagEig2','ePidx1','ePidx2','eP1','eP2'});

%% Initialisieren der Vektoren des for-loops

% Imaginaeranteil kontinuierlich steigend oder fallend
tmp = CharEx(:,6); % entspricht ImagEigSort(2) in BerechnungMathieu
tmpNeg = CharEx(:,5);  % entspricht ImagEigSort(1) in BerechnungMathieu

% Initialisieren der Ausgabevektoren
ImagEigCorrected = nan(size(tmp));
ImagEigCorrectedNeg = nan(size(tmp));

% Initialisieren des Buffers zum Erkennen des Erreichens von 90%
buffer = tmp(1); % Wird immer mit dem letzen Wert beschriebn
ImagEigCorrected(1) = tmp(1);

%% For-loop mit aus mat Datei geladenen, nicht korrigierten Werten
for idx = 2: length(tmp)
    if buffer <= tmp(idx) || (abs(tmp(idx)) < 10^-5) % Wert uebernehmen
        ImagEigCorrected(idx) = tmp(idx);  % steigender pos. Wert
        ImagEigCorrectedNeg(idx) = tmpNeg(idx); % fallender neg. Wert
        buffer = 0;
        bufferNeg = 0;
    else % 'korrigierten' Wert nehmen fuer kontinuierlichen Verlauf
        ImagEigCorrected(idx) = 2*buffer + tmpNeg(idx); % korrigierter pos. Wert
        ImagEigCorrectedNeg(idx) = 2*bufferNeg + tmp(idx); % korrigierter neg. Wert
    end
    % Buffer ueberschreiben mit aktuellem Wert
    buffer = max(tmp(idx),buffer); % Maximum speichern
    bufferNeg = min(tmpNeg(idx),bufferNeg); % Minimum speichern
end

%% test angle

               
%% Plot fuer mit atan berechnete nicht-korrigierte und korrigierte Werte
% Ueberpruefe die Verlaeufe mit atan nicht-korrigiert und korrigiert
cl = lines;
figure(1); 
set(gcf,'Name',"ImagEigCorrectedAtan")
subplot(2,1,1); plot(ImagEigCorrected,'-*',"Color",cl(1,:));
hold on; plot(tmp,'-*',"Color",cl(2,:));  grid on;
legend('ImagEigCorrected','Imag (atan)')
subplot(2,1,2); plot(ImagEigCorrectedNeg,'-*',"Color",cl(1,:)); 
hold on; plot(tmpNeg,'-*', "Color",cl(2,:));   grid on;

%% Berechnung der Winkel der char. Mulitplikatoren und der Quadranten

% Auslesen der charakteristische Multiplikatoren (Eigenwerte der Monodromiematrix)
T = 2*pi;
eP1 = tableCharEx.eP1; % in der von Matlab ausgegebenen Ordnung 
eP2 = tableCharEx.eP2;

phi = atan2(imag([eP1,eP2]),real([eP1,eP2])); % Winkelberechnung der char. Multiplikatoren  
phiDeg = phi *180/pi;
atan_eP = 1/T * phi;

% Berechnung der Quadranten
quadrant = 3 * ones(size(eP1)); % Default Wert 3. Quadrant
quadrant(real(eP1) >= 0 & imag(eP1) >= 0) = 1; % 1
quadrant(real(eP1) < 0 & imag(eP1) >= 0) = 2;
quadrant(real(eP1) >= 0 & imag(eP1) < 0) = 4;

quadrant2 = 3 * ones(size(eP1));
quadrant2(real(eP2) >= 0 & imag(eP2) >= 0) = 1;
quadrant2(real(eP2) < 0 & imag(eP2) >= 0) = 2;
quadrant2(real(eP2) >= 0 & imag(eP2) < 0) = 4;

atan_ePm(:,1) = mod(atan_eP(:,1),-0.5); % modulo
atan_ePm(:,2) = mod(atan_eP(:,2),0.5);

%%  Plot atan2 mit und ohne modulo und Quadranten
% Ueberpruefe die Verlaeufe mit atan2 mit modulo, negativ und positiv, und korrigierte
% Verlaeufe.

figure(2); set(gcf,'Name',"ImagEigAtan2Quadrant")
subplot(3,1,1); plot(atan_ePm(:,1),'-*',"Color",cl(1,:));
hold on; plot(atan_ePm(:,2),'-o',"Color",cl(2,:)); plot(ImagEigCorrectedNeg,'-x',"Color",cl(3,:));
legend('mod atan2 ePm1','mod atan2 ePm2','ImagEigCorrected')
subplot(3,1,2); plot(phiDeg(:,1),'-*');
hold on; plot(phiDeg(:,2),'-o');
legend('atan2 ePm1','atan2 ePm2')
subplot(3,1,3)
plot(quadrant,'*-'); hold on; plot(quadrant2,'*-');
legend('quadrant1','quadrant2')

%% Plot atan2 mit modulo
% Ueberpruefe die Verlaeufe mit atan2 mit modulo
figure(3); set(gcf,'Name',"ImagEigAtan2Modulo")
subplot(2,1,1); plot(mod(phiDeg(:,1),-180),'-*',"Color",cl(1,:)); 
hold on; plot(mod(phiDeg(:,2),180),'-o',"Color",cl(2,:));
subplot(2,1,2); plot(-mod(phiDeg(:,1),-180),'-*');

%% Plot Ergebnisse mit atan, atan2 und korrigiert
figure(4); 
set(gcf,'Name',"ImagEigCorrectedAtan2")

subplot(2,1,1);  plot(atan_ePm(:,2),'-o'); 
legend('Imag (atan2)')
subplot(2,1,2);  plot(tmp,'-o',"Color",cl(1,:)); hold on;  plot(ImagEigCorrected,'-*',"Color",cl(2,:));
legend('Imag (atan)','Corrected')

%% Nicht benutzter Code
% Auslesen der charakteristische Multiplikatoren (Eigenwerte der Monodromiematrix)
% ePidx1 = tableCharEx.ePidx1; % geordnet
% ePidx2 = tableCharEx.ePidx2;
% tan_ePidx1 = 1/T * atan(imag(ePidx1)./real(ePidx1)); % Berechnung mit atan
% tan_eP1 = 1/T * atan(imag(eP1)./real(eP1));
% tan_eP = 1/T * atan(imag([eP1,eP2])./real([eP1,eP2]));
% atan_ePidx1 = 1/T * atan2(imag(ePidx1),real(ePidx1)); % Berechnung mit atan2
% atan_eP1 = 1/T * atan2(imag(eP1),real(eP1));
% 
% figure(5); 
% subplot(4,1,1); plot(atan_ePm(:,1),'-*');
% hold on; plot(atan_ePm(:,2),'-o'); plot(ImagEigCorrectedNeg,'-x');
% legend('atan2 ePm1','atan2 ePm2','atan ePm1')
% hold on; plot(ImagEigCorrected,'-*'); 
% 
% subplot(4,1,2); plot(phiDeg(:,1),'-*');
% hold on; plot(phiDeg(:,2),'-o');
% subplot(4,1,3)
% plot(tan_eP(:,1),'-*');
% hold on; plot(tan_eP(:,2),'-o');
% subplot(4,1,4)
% plot(quadrant,'*-'); hold on; plot(quadrant2,'*-');

