%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Erzeugen des Stabilitaetsgebirges in den Grenzen von nu_02 und nu_C2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;

% Zeichne Netz mit ('EdgeColor',[0 0 0])
edgeColorOn = 'none'; %''; %;

% Zeichne Netz mit ('EdgeColor',[0 0 0])


% Schrittweite
SW = 0.1;

% Grenzen fuer nu_02 und nu_C2
unt0 = 0;
untC = 0;
ob0 = 20;
obC = 20;

% Parameter
D = 0.001;
t0 = 0.0;
T = 2*pi;

fDir = 'figureFolderStab'; % Ordner Abbildungen
if ~isdir(fDir) %#ok<ISDIR>
    mkdir(fDir)
end

dDir = 'dataFolder'; % Ordner mat-files
if ~isdir(dDir) %#ok<ISDIR>
    mkdir(dDir)
end

% Matlab mat file Name
matName = [strrep(sprintf('zMat_ob0%d_obC%d_SW%2.2f_D%2.1e', ob0, obC, SW, D),...
    '.','dot'),'.mat'];
fileName = fullfile(dDir,matName);
loadMat = 1; % Lade die Ergebnisse wenn sie bereits als mat Datei vorliegen

% Zaehlveriablen
oidx = 1;
AnzGl = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Vorbereitung der in den Schleifen zu fuellenden Arrays zwecks
% Programmbeschleunigung
Diagonal = diag(ones(AnzGl,1));
Monodromie = zeros(AnzGl);
CharEx = zeros(length(unt0:SW:ob0),AnzGl*3);
plotwert = zeros(length(unt0:SW:ob0)*length(untC:SW:obC),3);
Z = zeros(length(unt0:SW:ob0),length(untC:SW:obC));
Z0 = zeros(length(unt0:SW:ob0),length(untC:SW:obC));
Z1 = zeros(length(unt0:SW:ob0),length(untC:SW:obC));
Z2 = zeros(length(unt0:SW:ob0),length(untC:SW:obC));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[X,Y] = meshgrid(unt0:SW:ob0,untC:SW:obC);

if exist(fileName,'file') == 2 && loadMat == 1
    load(fileName,'Z','Z1','Z2','Z0'); % Laden der Z Matrizen
else
    for nu_02 = unt0:SW:ob0
        lidx = 1;
        for nu_C2 = untC:SW:obC

            %Optionen fuer die Genauigkeit und Toleranz fuer den ODE-Solver
            options = odeset('RelTol',1e-10,'AbsTol',1e-12);

            for k = 1:AnzGl
                sol = ode45(@(psi,x)MathieuDGL(psi,x,D,nu_02,nu_C2),[t0,T],Diagonal(:,k),options);
                MonoVek = deval(sol,T);
                Monodromie(:,k) = MonoVek;
            end

            % charakteristische Multiplikatoren (Eigenwerte der Monodromiematrix)
            p = [1, -(Monodromie(1,1) + Monodromie(2,2)), ...
                (Monodromie(1,1)*Monodromie(2,2) - Monodromie(1,2)*Monodromie(2,1))];
            eP = roots(p);
            realEP = real(eP); % Realanteil
            [~,idx] = sort(realEP); % Index Realanteil: Kleinerer Realanteil im 1. Element
            ePidx = eP(idx); % Geordneter Vektor

            ReLambda = 1/T * log(abs(eP(1))); % 1. Element ungeordnete Vektor
            ReLambda0 = 1/T * log(abs(eP(2))); % 2. Element ungeordnete Vektor

            ReLambda1 = 1/T * log(abs(ePidx(1))); % 1. Element geordneter Vektor
            ReLambda2 = 1/T * log(abs(ePidx(2))); % 2. Element geordneter Vektor

            Z(lidx,oidx) = ReLambda;
            Z0(lidx,oidx) = ReLambda0;
            Z1(lidx,oidx) = ReLambda1;
            Z2(lidx,oidx) = ReLambda2;
            lidx=lidx+1;
        end
        oidx = oidx+1;
    end
    save(fileName,'Z','Z1','Z2','Z0'); % Speichern der Z Matrizen
end

%% Grafische Darstellung

fs = 18;
figure(1), set(gcf,'Color',ones(1,3))
surf(X,Y,Z,'EdgeColor','none') %
view(30,35) %view([5 14])
xlim([0 ob0]);
ylim([0 obC]);
zlim([-2 2]);
xlabel('$\rm{Parameter} \; \nu_0^2 \rm{[-]}$','interpreter','latex','FontSize', fs);
ylabel('$\rm{Parameter} \; \nu_C^2 \rm{[-]}$','interpreter','latex','FontSize', fs);
zlabel('$Re(s_{R1}) \rm{[-]}$','interpreter','latex','FontSize', fs);

print(fullfile(fDir,'StabilitaetsgebirgePos_20_20'),'-dpng')


fs = 18;
figure(2)
surf(X,Y,Z0,'EdgeColor','none') %
view(-160,25) %view([5 14])
xlim([0 ob0]);
ylim([0 obC]);
zlim([-2 2]);
xlabel('$\rm{Parameter} \; \nu_0^2 \rm{[-]}$','interpreter','latex','FontSize', fs);
ylabel('$\rm{Parameter} \; \nu_C^2 \rm{[-]}$','interpreter','latex','FontSize', fs);
zlabel('$Re(s_{R2}) \rm{[-]}$','interpreter','latex','FontSize', fs);

print(fullfile(fDir,'StabilitaetsgebirgeNeg_20_20'),'-dpng')


%% Graphiken bis K0 = 2, KC = 2.5 und ohne negative Realteile
[Xlim2,Ylim2] = meshgrid(unt0:SW:2.5,untC:SW:2.5);
Zlim2 = Z(1:size(Xlim2,1),1:size(Xlim2,2));
Zlim20 = Z0(1:size(Xlim2,1),1:size(Xlim2,2));

figure(3)
surf(Xlim2,Ylim2,Zlim2)
view([27 24])
xlabel('$\rm{Parameter} \; \nu_0^2 \rm{[-]}$','interpreter','latex','FontSize', fs);
ylabel('$\rm{Parameter} \; \nu_C^2 \rm{[-]}$','interpreter','latex','FontSize', fs);
zlabel('$Re(s_{R1}) \rm{[-]}$','interpreter','latex','FontSize', fs);

print(fullfile(fDir,'StabilitaetsgebirgePos_2dot5'),'-dpng')


figure(4)
surf(Xlim2,Ylim2,Zlim20)
view(27,18)
xlabel('$\rm{Parameter} \; \nu_0^2 \rm{[-]}$','interpreter','latex','FontSize', fs);
ylabel('$\rm{Parameter} \; \nu_C^2 \rm{[-]}$','interpreter','latex','FontSize', fs);
zlabel('$Re(s_{R2}) \rm{[-]}$','interpreter','latex','FontSize', fs);

print(fullfile(fDir,'StabilitaetsgebirgeNeg_2dot5'),'-dpng')

%%

figure(5)
if isempty(edgeColorOn)
    surf(X,Y,Z1); view(30,35); fs2 = 2*fs; 
    set(gcf,"Position",1.0e+03 * [-2.5590    0.3710    2.5600    1.3273])
    set(gca,'FontSize',fs2)
else
  surf(X,Y,Z1,'EdgeColor','none','FaceAlpha', 0.8); view([27 24]); fs2 = fs;
  hold on;
  filtX = 1:5:length(X);
  surf(X(filtX,filtX),Y(filtX,filtX),Z1(filtX,filtX),'FaceColor','none');
end

xlim([0 ob0]);
ylim([0 obC]);
zlim([-2 2]);
xlabel('$\rm{Parameter} \; \nu_0^2 \rm{[-]}$','interpreter','latex','FontSize', fs2);
ylabel('$\rm{Parameter} \; \nu_C^2 \rm{[-]}$','interpreter','latex','FontSize', fs2);
zlabel('$Re(s_{max(real)}) \rm{[-]}$','interpreter','latex','FontSize', fs2);

print(fullfile(fDir,['StabilitaetsgebirgeMaxReal_20_20',edgeColorOn]),'-dpng')
 

figure(6)
if isempty(edgeColorOn)
    surf(X,Y,Z2); view(30,35); fs2 = 2*fs; 
    set(gcf,"Position",1.0e+03 * [-2.5590    0.3710    2.5600    1.3273])
    set(gca,'FontSize',fs2)
else
  %surf(X,Y,Z2,'EdgeColor','none'); view([5 14]);  fs2 = fs;
  surf(X,Y,Z2,'EdgeColor','none','FaceAlpha', 0.8); view([27 24]); fs2 = fs;
  hold on;
  filtX = 1:5:length(X);
  surf(X(filtX,filtX),Y(filtX,filtX),Z2(filtX,filtX),'FaceColor','none');

end

xlim([0 ob0]);
ylim([0 ob0]);
zlim([-2 2]);
xlabel('$\rm{Parameter} \; \nu_0^2 \rm{[-]}$','interpreter','latex','FontSize', fs2);
ylabel('$\rm{Parameter} \; \nu_C^2 \rm{[-]}$','interpreter','latex','FontSize', fs2);
zlabel('$Re(s_{min(real)}) \rm{[-]}$','interpreter','latex','FontSize', fs2);

print(fullfile(fDir,['StabilitaetsgebirgeMinReal_20_20',edgeColorOn]),'-dpng')
 


