%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Berechnung der Struttschen Karten in den Grenzen von nu_02 und nu_C2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

loadMat = 0;  % mat-file laden, wenn Ergebnisse mit gleichem D und tSpan vorhanden
useTspan = 0; % tspan statt [t0,T] in der ode45 benutzt

%Schrittweite
SW = 0.1;

% Grenzen fuer nu_02 und nu_C2
unt0 = 0;
untC = 0;
ob0 = 9;
obC = 9;

% Zaehlvariablen
lidx = 1; % Index charakteristische Multiplikatoren
oidx = 1; % Index charakteristische Exponenten

% Anzahl der Gleichungen des DGL-Systems
Nz = 2;

% Parameter
D = 0.15; %0.001; % 0.3 0.5
t0 = 0.0;
T = 2*pi;
Dt = 0.0001;
tspan = t0:Dt:(T+Dt);

%Anfangsbedingungen
Diagonal = diag(ones(Nz,1));


% Matlab mat-file Name
matName = [strrep(sprintf('STRUTTscheKarte_D%2.1e_tSpanUsed%d',D,useTspan),'.','dot'),...
    '.mat'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Vorbereitung der in den Schleifen zu fuellenden Arrays zwecks
% Programmbeschleunigung
lenNu02 = length(unt0:SW:ob0); % Anzahl Nu_02 Werte
lenNuC2 = length(untC:SW:obC); % Anzahl Nu_C2 Werte
lenNu = lenNu02 * lenNuC2; % Anzahl Kombination Nu_C2 Werte
lenNuDiag = min(lenNu02,lenNuC2); % Anzahl Werte nu_02 == nu_C2
Monodromie = zeros(Nz);
CharMult =  zeros(lenNu,2);
CharEx = zeros(lenNuDiag,Nz*4+2+2); % nu, nc, Real1,2, Imag12, Imag1,2, Pole
plotwertstabil = zeros(lenNu,3);
plotwertinstabil = zeros(lenNu,3);
compareMonoVek = zeros(2 + Nz*2, lenNuDiag);

%nuCSwitchVec = [0.125,1,1.5^2,2^2,2.5^2] - 0.1;
nuCSwitchVec = [0.2,1,1.5^2,2^2,2.5^2] - 0.15;
n = 0;
cntN = 1;

buffer = 0;
bufferNeg = 0;
vecSwitch = zeros(lenNu02,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist(matName,'file') == 2 && loadMat == 1
    load(matName,'Char*','plotwert*', 'nu*','compare*');

else

    for nu_02 = unt0:SW:ob0

        for nu_C2  = untC:SW:obC

            %Schleife intergiert Nz Mal numerisch nacheinaner fuer die
            %Spaltenvektren der Einheitsmatrix, Auswertung der Loesung bei T
            %wird in MonoVek geschrieben und diese zur Monodromiematrix
            %zusammengesetzt
            options = odeset('RelTol',1e-10,'AbsTol',1e-12);

            for k = 1 : Nz
                sol = ode45(@(psi,x)MathieuDGL(psi,x,D,nu_02,nu_C2),[t0,T],Diagonal(:,k),options);
                MonoVek  = deval(sol,T);

                if useTspan == 1
                    sol1 = ode45(@(psi,x)MathieuDGL(psi,x,D,nu_02,nu_C2),tspan,Diagonal(:,k),options); %
                    MonoVek1 = deval(sol1,T);

                    if any( abs(MonoVek - MonoVek1)./abs(MonoVek) > 10^(-2) ) % Differenz groesser als 1%
                        disp('MonoVek not equal')
                        compareMonoVek(:,lidx) = [nu_02; nu_C2; MonoVek;MonoVek1];
                    end

                    Monodromie(:,k) = MonoVek1;
                else
                    Monodromie(:,k) = MonoVek;
                end
            end

            % charakteristische Multiplikatoren (Eigenwerte der Monodromiematrix)
            p = [1, -(Monodromie(1,1) + Monodromie(2,2)), ...
                Monodromie(1,1)*Monodromie(2,2)-Monodromie(1,2)*Monodromie(2,1)];
            eP = roots(p)';
            realEP = real(eP); % Realanteil
            [~,idx] = sort(realEP); % Index des Realanteil: Kleinerer Realanteil im 1. Element
            ePidx = eP(idx);
            CharMult(lidx,:) = ePidx; % Spaltenvektor von Matlab als Zeilenvektor einsortiert

            % Betrag der Eigenwerte zur Stabilitaetsanalyse
            eBetr = abs(ePidx);

            % Charakteristische Exponenten
            if nu_C2 == nu_02
                RealEig = 1/T * log(eBetr);
                ImagEig = 1/T * atan(imag(ePidx)./real(ePidx));
                [~,idx] = sort(ImagEig);
                ImagEigSort = ImagEig(idx);


                % Imaginaeranteil kontinuierlich steigend oder fallend
                tmp = ImagEigSort(2);
                tmpNeg = ImagEigSort(1);

                if buffer <= tmp || (abs(tmp) < 10^-5) % Wert uebernehmen
                    ImagEigCorrected = tmp; % steigender pos. Wert
                    ImagEigCorrectedNeg = tmpNeg;  % fallender neg. Wert
                    buffer = 0;
                    bufferNeg = 0;
                else  % 'korrigierten' Wert nehmen fuer kontinuierlichen Verlauf
                    % Wird der 'else'-Zweig getriggert, ist der buffer
                    % bei 90 Grad (1.57 rad). Der positive Imaginaerteil ist damit
                    % die Summe aus 180 Grad und dem negativn Winkel, dessen
                    % Wert von -90 Grad zu 0 Grad laeuft.
                    ImagEigCorrected = 2*buffer + tmpNeg; % korrigierter pos. Wert
                    ImagEigCorrectedNeg = 2*bufferNeg + tmp; % korrigierter neg. Wert
                    vecSwitch(oidx)= 1;

                end
                % Buffer ueberschreiben mit aktuellem Wert
                buffer = max(tmp,buffer); % Maximum speichern
                bufferNeg = min(tmpNeg,bufferNeg); % Minimum speichern

                % Additionsterm n*2*pi/T
                if cntN <= length(nuCSwitchVec) && ... % Sicherheitscheck, damit n nicht groesser wird als die Anzahl der 'Switchstellen'
                        nu_C2 > nuCSwitchVec(cntN) && ... % n nur bei erreichen der Switchstellen umstellen
                        abs(RealEig(1) - RealEig(2))> eps % die Realteile muessen gleich sein
                    n =  n + 0.5;
                    cntN = cntN+1;
                end

                nAdd = n*2*pi/T;
                ImagEigSortN = [ImagEigCorrectedNeg, ImagEigCorrected] + [-nAdd,nAdd];


                CharEx(oidx,:) = [nu_02, nu_C2, RealEig, ...
                    min(ImagEigSort), max(ImagEigSort), ImagEigSortN, ePidx, eP];

                oidx = oidx + 1;

            end

            % 1=stabil , -1=instabil Anmerkung nicht vergessen (zu =1)
            % sortieren der Parameterkombinationen in Array mit stabilen und
            % mit instabilen Kombinationen von nu_02 und nu_C2 anhand des
            % charakteristischen Mutiplikators
            % b = 0;

            if eBetr(1)<1 && eBetr(2)<1
                b = 1;
                plotwertstabil(lidx,:) = [nu_02,nu_C2,b];
            elseif eBetr(1)>1 || eBetr(2)>1
                b = -1;
                plotwertinstabil(lidx,:) = [nu_02,nu_C2,b];
            end
            lidx = lidx + 1;

        end
    end

    %Loeschen der reinen 0-Zeilen
    plotwertstabil = plotwertstabil(any(plotwertstabil,2),:);
    plotwertinstabil = plotwertinstabil(any(plotwertinstabil,2),:);

    save(matName,'Char*','plotwert*', 'nu*','compare*');
end

tableCharEx = array2table(CharEx,'VariableNames',...
    {'nu_02','nu_C2','RealEig1','RealEig2','ImagEig1_0','ImagEig2_0',...
    'ImagEig1','ImagEig2','ePidx1','ePidx2','eP1','eP2'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Grafische Darstellung: Struttsche Karte  (Plotten aller Wertepaare
%nu_02/nu_C2, die zu stabiler Loesung fuehren) sowie Realteile der
% charakteristischen Exponenten
cl = lines;
hf = figure(2);
hf.Position = [10 10 600 600];
tiledlayout(11,1);

h(1) = nexttile([7 1]);
% Struttsche Karte
scatter(plotwertstabil(:,1),plotwertstabil(:,2),1)
ylabel('$\rm{Parameter} \; \nu_C^2 \;\;\rm{[-]}$','interpreter','latex','Position', [-0.5 nu_C2/2],'FontSize', 10);
grid on

% Realteile char. Exponenten
h(2) = nexttile([2 1]);

xachse = unt0:SW:ob0;
plot(xachse,CharEx(:,3:4))

ylabel('$Re(s_R) \;\; \rm{[-]}$','interpreter','latex','Position', [-0.5 -D],'FontSize', 10);
grid on;

tmp = tableCharEx.ImagEig2;
tmp1 = nan(size(tmp));
buffer = tableCharEx.ImagEig2(1);
tmp1(1) = tableCharEx.ImagEig2(1);
for idx = 2: length(tmp)
    tmp1(idx) = max(buffer,tmp(idx));
    buffer = tmp1(idx);
end
CharEx(:,8) = tmp1;

tmp = tableCharEx.ImagEig1;
tmp1 = nan(size(tmp));
buffer = tableCharEx.ImagEig1(1);
tmp1(1) = tableCharEx.ImagEig1(1);
for idx = 2: length(tmp)
    tmp1(idx) = min(buffer,tmp(idx));
    buffer = tmp1(idx);
end
CharEx(:,7) = tmp1;

h(3) = nexttile([2 1]);
plot(xachse,CharEx(:,7:8)); hold on;
% plot(xachse,CharEx(:,5),'--', 'color',cl(1,:))
% plot(xachse,CharEx(:,6),'--', 'color',cl(2,:))
grid on;
ylim(gca,[-0.15,3]);

% h(4) = nexttile([2 1]);
% plot(xachse,real(tableCharEx.ePidx1),xachse,real(tableCharEx.ePidx2)); hold on;
%
% h(5) = nexttile([2 1]);
% plot(xachse,imag(tableCharEx.ePidx1),xachse,imag(tableCharEx.ePidx2)); hold on;

grid on;

xlabel('$\rm{Parameter} \; \nu_0^2 \;\;\rm{[-]}$','interpreter','latex','FontSize', 10);
ylabel('$Im(s_R) \;\; \rm{[-]}$','interpreter','latex','Position', [-0.5 -D],'FontSize', 10);
linkaxes(h,'x')

pngname = strrep(matName,'.mat','');
print(pngname, '-dpng')

% figure;
% for idx = 1: 6
%     vec1 = (1+(idx-1)*13):13*idx; plot(real(tableCharEx.ePidx1(vec1)), imag(tableCharEx.ePidx1(vec1)),'o'); hold on;
% end

return;
%% Zusaetzliche Graphiken

figure(2)
% Darstellung des Real- und Imaginaeranteile der Pole
fs = 16;
CharExReal = real(CharEx(:,7:8));
CharExImag = imag(CharEx(:,7:8));
tiledlayout(2,1);
nexttile
plot(xachse,CharExReal)
ylabel('Realanteil Pole','interpreter','latex','Position', [-0.5 -D],'FontSize',fs);

nexttile
plot(xachse,CharExImag)
xlabel('$\rm{Parameter} \; \nu_0^2 \;\;\rm{[-]}$','interpreter','latex','FontSize', fs);
ylabel('Imaginaeranteil Pole','interpreter','latex','Position', [-0.5 -D],'FontSize', fs);

figure(3);
% 'Inverse' Struttsche Karte
scatter(plotwertinstabil(:,1),plotwertinstabil(:,2),1)
ylabel('$\rm{Parameter} \; \nu_C^2 \;\;\rm{[-]}$','interpreter','latex','Position', [-0.5 nu_C2/2],'FontSize', 10);


%% Code zum Testen (liefert Abweichung)
% [t45,y45]= ode45(@(psi,x)MathieuDGL(psi,x,D,nu_02,nu_C2),tspan,Diagonal(:,k));
%  MonoVek2 = y45(abs(t45 - T) < Dt/2,:)';