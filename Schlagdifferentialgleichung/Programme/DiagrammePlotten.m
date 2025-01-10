%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Der Code in dieser Datei ist auf die bearbeiteten Faelle zugeschnitten, um
% diese fuer die Arbeit aufzubereiten. Sie funktionieren teils nur fuer
% bestimmte Parameterbereiche. Um allerdings die Diagramme im Schwebeflug,
% die Eigenwerte in der komplexen Ebene und die Real- und Imaginaerteile
% anzuzeigen, reicht der Code des ersten Abschnittes. Um einen einzelnen
% Abschnitt auszufuehren muss sich der Cursor in diesem Abschnitt befinden
% und Strg+Enter gedrueckt werden.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% allgemein gueltiger Code zum Plotten aller Faelle
clc; clear; close all

load("Workspace.mat","damp","freq","MuMin","SW","MuMax", ...
    "CharExRe","CharExIm","nu0","Blatt")

xachse = MuMin:SW:MuMax; % Diagrammgrenzen
nFig = 1; % Nummer erste Figure
useCorVal = 1; % Ueberschreiben der Realteil


%% Glaette den Verlauf des Realteils
CharExRe1 = CharExRe(:,1);
CharExRe2 = CharExRe(:,2);

CharExRePos = max(CharExRe1,CharExRe2); % Positiver Eigenwert: Glatter Verlauf
CharExRe1_NegIdx =  abs(CharExRePos - CharExRe1) > eps; % Index: Negativer Wert 1. Realteil 
offset = CharExRe1(1); % Offset 
CharExReNeg =  - (CharExRePos - offset) + offset; % Negativer EW: Gespiegelter positiver EW

CharExRe1Cor = CharExRe1; % Korrigierter 1. Eigenwert
CharExRe1Cor(CharExRe1_NegIdx) = CharExReNeg(CharExRe1_NegIdx); % Ersetze negative Werte durch gespiegelte positive EW

CharExRe2Cor = CharExRe2; % Korrigierter 1. Eigenwert
CharExRe2Cor(~CharExRe1_NegIdx) = CharExReNeg(~CharExRe1_NegIdx); % Ersetze negative Werte durch gespiegelte positive EW

if useCorVal == 1

    CharExRe(:,1) = CharExRe1Cor;
    CharExRe(:,2) = CharExRe2Cor;
    CharExIm(:,2) = - CharExIm(:,1);
end


%% Plot Real- und Imaginaerteils
figure(nFig); nFig = nFig + 1;

ax1(1) = subplot(2,1,1);
plot(xachse,CharExRe,'LineWidth',1.5,'Color','k');
% set(gca,'ydir', 'reverse');
grid on;
xlabel('$\mu$','interpreter','latex','FontSize', 14);
ylabel('${Re(s_R)}$','interpreter','latex','FontSize', 14);

ax1(2) = subplot(2,1,2);
plot(xachse,CharExIm,'LineWidth',1.5,'Color','k');
grid on;
ylim('padded')
xlabel('$\mu$','interpreter','latex','FontSize', 14);
ylabel('${Im(s_R)}$','interpreter','latex','FontSize', 14);
linkaxes(ax1,'x')

saveas(gcf,[pwd '/Plots/RealImaginaerteile.png'])
saveas(gcf,[pwd '/Plots/RealImaginaerteile.fig'])


%% Plot Zustaende Schwebeflug (damping und frequency)
figure(nFig); nFig = nFig + 1;
set(gcf,'Position', [100 100 700 700]); 
s = scatter(damp,freq,100,'linewidth',2);
m = s.Marker;
s.Marker = 'x';
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
set(gca,'XLim',[-2*nu0 2*nu0])
xlabel('${Re(s_R)}$','interpreter','latex','FontSize', 14);
ylabel('${Im(s_R)}$','interpreter','latex','FontSize', 14);
saveas(gcf,[pwd '/Plots/Schwebeflug.fig'])
saveas(gcf,[pwd '/Plots/Schwebeflug.png'])


%% Plotten der Eigenwerte in der komplexen Ebene
figure;
scatter(CharExRe,CharExIm,1);

% for idx = [1,30:10:50,100:100:700] %length(CharExRe)
%     text(CharExRe(idx), CharExIm(idx),num2str(idx));
% end

saveas(gcf,[pwd '/Plots/Eigenwerte.fig'])
saveas(gcf,[pwd '/Plots/Eigenwerte.png'])

if Blatt == 1
    return;
end


figure;
scatter(CharExRe(:,1),CharExIm(:,1),'.');


% for idx = [1,30:10:50,100:100:700] %length(CharExRe)
% 
%     text(CharExRe(idx,1), CharExIm(idx,1),num2str(idx,1));
% end
% 
% figure;
% scatter(CharExRe(:,2),CharExIm(:,2),'.');
% 
% 
% for idx = [1,30:10:50,100:100:700] %length(CharExRe)
% 
%     text(CharExRe(idx), CharExIm(idx),num2str(idx));
% end


%% 3-Blatt-Rotor Zentrales Schlaggelenk gamma1 = 5,858

% load("Workspace.mat")
% close all

figure(nFig); nFig = nFig + 1;
set(gcf,'Position', [100 100 700 700]);
s = scatter(damp,freq,100,'linewidth',2);
m = s.Marker;
s.Marker = 'x';
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
set(gca,'XLim',[-2*nu0 2.5*nu0],'YLim',[-2.1*nu0 2.1*nu0])
% viscircles([0 0],nu0,'color','k','linewidth',0.8)
line([0 0.75],[0 sqrt(nu0^2-0.75^2)])
xlabel('${Re(s_R)}$','interpreter','latex','FontSize', 14);
ylabel('${Im(s_R)}$','interpreter','latex','FontSize', 14);
daspect([1 1 1])
saveas(gcf,[pwd '/Plots/Schwebeflug.fig'])
saveas(gcf,[pwd '/Plots/Schwebeflug.png'])


%Diagrammgrenzen
xunten = min(CharExRe(length(CharExRe),AnzGl/2))-0.2;
xoben = max(CharExRe(length(CharExRe),AnzGl))+0.2;

xachse = MuMin:SW:MuMax;


% Plotten des Realteils
figure;
plot(xachse,CharExRe,'LineWidth',1.5,'Color','k');
set(gca,'ydir', 'reverse');
xlabel('$\mu$','interpreter','latex','FontSize', 14);
ylabel('${Re(s_R)}$','interpreter','latex','FontSize', 14);
saveas(gcf,[pwd '/Plots/Realteile.fig'])
saveas(gcf,[pwd '/Plots/Realteile.png'])

if size(CharExIm,2) <5
    return;
end


% Plotten des Imaginaerteils
figure('Position', [10 10 600 300]);
plot(xachse,CharExIm(:,5),'LineWidth',1.5,'Color','k');
ylim('padded')
xlabel('$\mu$','interpreter','latex','FontSize', 14);
ylabel('${Im(s_R)}$','interpreter','latex','FontSize', 14);
saveas(gcf,[pwd '/Plots/Imaginaerteile.fig'])
saveas(gcf,[pwd '/Plots/Imaginaerteile.png'])


PlotRe1 = [CharExRe(:,4),CharExRe(:,3)];
PlotIm1 = [CharExIm(:,4),CharExIm(:,3)];


figure('Position', [10 10 600 AnzGl*2*50]);
tl = tiledlayout(AnzGl*2,1);

nexttile([2 1]);
scatter(CharExRe(:,6),CharExIm(:,6),1,'k')
set(gca,'XTick',[],'YTick',[1.5 round(freq(1),2) 2],'XColor','none','XLim',[xunten xoben])
ylim('padded')

nexttile([2 1]);
scatter(CharExRe(:,5),CharExIm(:,5),1,'k')
set(gca,'XTick',[],'YTick',[round(freq(2),2) 1],'XColor','none','XLim',[xunten xoben])
ylim('padded')

nexttile([4 1]);
scatter(PlotRe1,PlotIm1,1,'k')
set(gca,'XTick',[],'YTick',[round(freq(4),2) 0 round(freq(3),2)],'XColor','none','XLim',[xunten xoben])
ylim('padded')

nexttile([2 1]);
scatter(CharExRe(:,2),CharExIm(:,2),1,'k')
set(gca,'XLim',[xunten xoben],'YTick',[-1 round(freq(5),2)],'XColor','none')
ylim('padded')

nexttile([2 1]);
scatter(CharExRe(:,1),CharExIm(:,1),1,'k')
set(gca,'YTick',[-2 round(freq(6),2)],'XLim',[xunten xoben])
ylim('padded')

xlabel(tl,'${Re(s_R)}$','interpreter','latex','FontSize', 14);
ylabel(tl,'${Im(s_R)}$','interpreter','latex','FontSize', 14);

saveas(gcf,[pwd '/Plots/Eigenwerte.fig'])
saveas(gcf,[pwd '/Plots/Eigenwerte.png'])


%% 3-Blatt-Rotor Zentrales Schlaggelenk gamma2 = 12

load("Workspace.mat")
close all

figure('Position', [100 100 700 700]);
s = scatter(damp,freq,100,'linewidth',2);
m = s.Marker;
s.Marker = 'x';
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
set(gca,'XLim',[-2*nu0 2.5*nu0],'YLim',[-2.1*nu0 2.1*nu0])
viscircles([0 0],nu0,'color','k','linewidth',0.8)
line([0 0.75],[0 sqrt(nu0^2-0.75^2)])
xlabel('${Re(s_R)}$','interpreter','latex','FontSize', 14);
ylabel('${Im(s_R)}$','interpreter','latex','FontSize', 14);
daspect([1 1 1])
saveas(gcf,[pwd '/Plots/Schwebeflug.fig'])
saveas(gcf,[pwd '/Plots/Schwebeflug.png'])

xachse = MuMin:SW:MuMax;

%Diagrammgrenzen
xunten = min(CharExRe(length(CharExRe),AnzGl/2))-0.2;
xoben = max(CharExRe(length(CharExRe),AnzGl))+0.2;

% Plotten des Realteils
figure;
plot(xachse,CharExRe,'LineWidth',1.5,'Color','k');
set(gca,'ydir', 'reverse');
xlabel('$\mu$','interpreter','latex','FontSize', 14);
ylabel('${Re(s_R)}$','interpreter','latex','FontSize', 14);
saveas(gcf,[pwd '/Plots/Realteile.fig'])
saveas(gcf,[pwd '/Plots/Realteile.png'])


% Plotten des Imaginaerteils
figure('Position', [10 10 600 300]);
plot(xachse,CharExIm(:,5),'LineWidth',1.5,'Color','k');
ylim('padded')
xlabel('$\mu$','interpreter','latex','FontSize', 14);
ylabel('${Im(s_R)}$','interpreter','latex','FontSize', 14);
saveas(gcf,[pwd '/Plots/Imaginaerteile.fig'])
saveas(gcf,[pwd '/Plots/Imaginaerteile.png'])


PlotRe1 = [CharExRe(:,4:5),CharExRe(:,2:3)];
PlotIm1 = [CharExIm(:,4:5),CharExIm(:,2:3)];

figure('Position', [10 10 600 AnzGl*2*50]);
tl = tiledlayout(AnzGl*2,1);

nexttile([2 1]);
scatter(CharExRe(:,6),CharExIm(:,6),1,'k')
set(gca,'XTick',[],'YTick',[1.5 round(freq(1),2) round(max(CharExIm(:,6)),2)],'XColor','none','XLim',[xunten xoben])
ylim('padded')

nexttile([8 1]);
scatter(PlotRe1,PlotIm1,1,'k')
set(gca,'XTick',[],'YTick',[-1 round(freq(5),2) -0.5 round(freq(4),2) 0 round(freq(3),2) 0.5 round(freq(2),2) 1],'XColor','none','XLim',[xunten xoben])
ylim('padded')


nexttile([2 1]);
scatter(CharExRe(:,1),CharExIm(:,1),1,'k')
set(gca,'XLim',[xunten xoben],'YTick',[round(min(CharExIm(:,1)),2) round(freq(6),2) -1.5])
ylim('padded')

xlabel(tl,'${Re(s_R)}$','interpreter','latex','FontSize', 14);
ylabel(tl,'${Im(s_R)}$','interpreter','latex','FontSize', 14);


saveas(gcf,[pwd '/Plots/Eigenwerte.fig'])
saveas(gcf,[pwd '/Plots/Eigenwerte.png'])



%% 3-Blatt-Rotor Voll-Gelenkig

load("Workspace.mat")
close all

figure('Position', [100 100 700 700]);
s = scatter(damp,freq,100,'linewidth',2);
m = s.Marker;
s.Marker = 'x';
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
set(gca,'XLim',[-2*nu0 2.5*nu0],'YLim',[-2.1*nu0 2.1*nu0])
viscircles([0 0],nu0,'color','k','linewidth',0.8)
line([0 0.75],[0 sqrt(nu0^2-0.75^2)])
xlabel('${Re(s_R)}$','interpreter','latex','FontSize', 14);
ylabel('${Im(s_R)}$','interpreter','latex','FontSize', 14);
daspect([1 1 1])
saveas(gcf,[pwd '/Plots/Schwebeflug.fig'])
saveas(gcf,[pwd '/Plots/Schwebeflug.png'])

%Diagrammgrenzen
xunten = min(CharExRe(length(CharExRe),AnzGl/2))-0.2;
xoben = max(CharExRe(length(CharExRe),AnzGl))+0.2;

xachse = MuMin:SW:MuMax;


% Plotten des Realteils
figure;
plot(xachse,CharExRe,'LineWidth',1.5,'Color','k');
set(gca,'ydir', 'reverse');
xlabel('$\mu$','interpreter','latex','FontSize', 14);
ylabel('${Re(s_R)}$','interpreter','latex','FontSize', 14);
saveas(gcf,[pwd '/Plots/Realteile.fig'])
saveas(gcf,[pwd '/Plots/Realteile.png'])


% Plotten des Imaginaerteils
figure('Position', [10 10 600 300]);
plot(xachse,CharExIm(:,5),'LineWidth',1.5,'Color','k');
ylim('padded')
xlabel('$\mu$','interpreter','latex','FontSize', 14);
ylabel('${Im(s_R)}$','interpreter','latex','FontSize', 14);
saveas(gcf,[pwd '/Plots/Imaginaerteile.fig'])
saveas(gcf,[pwd '/Plots/Imaginaerteile.png'])


PlotRe1 = [CharExRe(:,4),CharExRe(:,3)];
PlotIm1 = [CharExIm(:,4),CharExIm(:,3)];


figure('Position', [10 10 600 AnzGl*2*50]);
tl = tiledlayout(AnzGl*2,1);

nexttile([2 1]);
scatter(CharExRe(:,6),CharExIm(:,6),1,'k')
set(gca,'XTick',[],'YTick',[1.5 round(freq(1),2) 2],'XColor','none','XLim',[xunten xoben])
ylim('padded')

nexttile([2 1]);
scatter(CharExRe(:,5),CharExIm(:,5),1,'k')
set(gca,'XTick',[],'YTick',[round(freq(2),2) 1],'XColor','none','XLim',[xunten xoben])
ylim('padded')

nexttile([4 1]);
scatter(PlotRe1,PlotIm1,1,'k')
set(gca,'XTick',[],'YTick',[round(freq(4),2) 0 round(freq(3),2)],'XColor','none','XLim',[xunten xoben])
ylim('padded')

nexttile([2 1]);
scatter(CharExRe(:,2),CharExIm(:,2),1,'k')
set(gca,'XLim',[xunten xoben],'YTick',[-1 round(freq(5),2)],'XColor','none')
ylim('padded')

nexttile([2 1]);
scatter(CharExRe(:,1),CharExIm(:,1),1,'k')
set(gca,'YTick',[-2 round(freq(6),2)],'XLim',[xunten xoben])
ylim('padded')

xlabel(tl,'${Re(s_R)}$','interpreter','latex','FontSize', 14);
ylabel(tl,'${Im(s_R)}$','interpreter','latex','FontSize', 14);


saveas(gcf,[pwd '/Plots/Eigenwerte.fig'])
saveas(gcf,[pwd '/Plots/Eigenwerte.png'])

%% 4-Blatt-Rotor Voll-Gelenkig

load("Workspace.mat")
close all

figure('Position', [100 100 700 700]);
s = scatter(damp,freq,100,'linewidth',2);
m = s.Marker;
s.Marker = 'x';
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
set(gca,'XLim',[-2*nu0 2.5*nu0],'YLim',[-2.1*nu0 2.1*nu0])
viscircles([0 0],nu0,'color','k','linewidth',0.8)
line([0 0.75],[0 sqrt(nu0^2-0.75^2)])
xlabel('${Re(s_R)}$','interpreter','latex','FontSize', 14);
ylabel('${Im(s_R)}$','interpreter','latex','FontSize', 14);
daspect([1 1 1])
saveas(gcf,[pwd '/Plots/Schwebeflug.fig'])
saveas(gcf,[pwd '/Plots/Schwebeflug.png'])

%Diagrammgrenzen
xunten = min(CharExRe(length(CharExRe),AnzGl/2))-0.2;
xoben = max(CharExRe(length(CharExRe),AnzGl))+0.2;

xachse = MuMin:SW:MuMax;

% Plotten des Realteils
figure;
plot(xachse,CharExRe,'LineWidth',1.5,'Color','k');
set(gca,'ydir', 'reverse','YLim',[xunten xoben]);
xlabel('$\mu$','interpreter','latex','FontSize', 14);
ylabel('${Re(s_R)}$','interpreter','latex','FontSize', 14);
saveas(gcf,[pwd '/Plots/Realteile.fig'])
saveas(gcf,[pwd '/Plots/Realteile.png'])


% Plotten des Imaginaerteils
figure('Position', [10 10 600 300]);
plot(xachse,CharExIm(:,6),'LineWidth',1.5,'Color','k');
%scatter(xachse,CharExIm,1)
ylim('padded')
xlabel('$\mu$','interpreter','latex','FontSize', 14);
ylabel('${Im(s_R)}$','interpreter','latex','FontSize', 14);
saveas(gcf,[pwd '/Plots/Imaginaerteile.fig'])
saveas(gcf,[pwd '/Plots/Imaginaerteile.png'])


PlotRe1 = [CharExRe(:,4),CharExRe(:,5)];
PlotIm1 = [CharExIm(:,4),CharExIm(:,5)];


figure('Position', [10 10 600 6*2*50]);
tl = tiledlayout(6*2,1);

nexttile([2 1]);
scatter(CharExRe(:,8),CharExIm(:,8),1,'k')
set(gca,'XTick',[],'YTick',[round(freq(1),2) 2],'XColor','none','XLim',[xunten xoben])
ylim('padded')

nexttile([2 1]);
scatter(CharExRe(:,7),CharExIm(:,7),1,'k')
set(gca,'XTick',[],'YTick',[round(freq(2),2) 1],'XColor','none','XLim',[xunten xoben])
ylim('padded')

nexttile([4 1]);
scatter(PlotRe1,PlotIm1,1,'k')
set(gca,'XTick',[],'YTick',[round(freq(5),2) 0 round(freq(4),2)],'XColor','none','XLim',[xunten xoben])
ylim('padded')

nexttile([2 1]);
scatter(CharExRe(:,2),CharExIm(:,2),1,'k')
set(gca,'XLim',[xunten xoben],'YTick',[-1 round(freq(7),2)],'XColor','none')
ylim('padded')

nexttile([2 1]);
scatter(CharExRe(:,1),CharExIm(:,1),1,'k')
set(gca,'YTick',[-2 round(freq(8),2)],'XLim',[xunten xoben])
ylim('padded')


xlabel(tl,'${Re(s_R)}$','interpreter','latex','FontSize', 14);
ylabel(tl,'${Im(s_R)}$','interpreter','latex','FontSize', 14);


saveas(gcf,[pwd '/Plots/Eigenwerte.fig'])
saveas(gcf,[pwd '/Plots/Eigenwerte.png'])


%% 5-Blatt-Rotor Voll-Gelenkig

load("Workspace.mat")
close all


figure('Position', [100 100 700 700]);
s = scatter(damp,freq,100,'linewidth',2);
m = s.Marker;
s.Marker = 'x';
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
set(gca,'XLim',[-2*nu0 2.5*nu0],'YLim',[-3.1*nu0 3.1*nu0])
viscircles([0 0],nu0,'color','k','linewidth',0.8)
line([0 0.75],[0 sqrt(nu0^2-0.75^2)])
xlabel('${Re(s_R)}$','interpreter','latex','FontSize', 14);
ylabel('${Im(s_R)}$','interpreter','latex','FontSize', 14);
daspect([1 1 1])
saveas(gcf,[pwd '/Plots/Schwebeflug.fig'])
saveas(gcf,[pwd '/Plots/Schwebeflug.png'])

%Diagrammgrenzen
xunten = min(CharExRe(:,1))-0.2;
xoben = max(CharExRe(:,AnzGl))+0.2;

xachse = MuMin:SW:MuMax;


% Plotten des Realteils
figure;
plot(xachse,CharExRe,'LineWidth',1.5,'Color','k');
set(gca,'ydir', 'reverse','YLim',[xunten xoben]);
xlabel('$\mu$','interpreter','latex','FontSize', 14);
ylabel('${Re(s_R)}$','interpreter','latex','FontSize', 14);
saveas(gcf,[pwd '/Plots/Realteile.fig'])
saveas(gcf,[pwd '/Plots/Realteile.png'])


% Plotten des Imaginaerteils
figure('Position', [10 10 600 300]);
plot(xachse,CharExIm(:,8),'LineWidth',1.5,'Color','k');
%scatter(xachse,CharExIm,1)
ylim('padded')
xlabel('$\mu$','interpreter','latex','FontSize', 14);
ylabel('${Im(s_R)}$','interpreter','latex','FontSize', 14);
saveas(gcf,[pwd '/Plots/Imaginaerteile.fig'])
saveas(gcf,[pwd '/Plots/Imaginaerteile.png'])


PlotRe1 = [CharExRe(:,5),CharExRe(:,8)];
PlotIm1 = [CharExIm(:,5),CharExIm(:,8)];
PlotRe2 = [CharExRe(:,4),CharExRe(:,7)];
PlotIm2 = [CharExIm(:,4),CharExIm(:,7)];
PlotRe3 = [CharExRe(:,3),CharExRe(:,6)];
PlotIm3 = [CharExIm(:,3),CharExIm(:,6)];


figure('Position', [10 10 600 10*2*50]);
tl = tiledlayout(10*2,1);

nexttile([2 1]);
scatter(CharExRe(:,10),CharExIm(:,10),0.8,'k')
set(gca,'XTick',[],'YTick',[round(freq(1),2) 3],'XColor','none','XLim',[xunten xoben])
ylim('padded')

nexttile([2 1]);
scatter(CharExRe(:,9),CharExIm(:,9),0.8,'k')
set(gca,'XTick',[],'YTick',[round(freq(2),2) 2],'XColor','none','XLim',[xunten xoben])
ylim('padded')

nexttile([4 1]);
scatter(PlotRe1,PlotIm1,1,'k')
set(gca,'XTick',[],'YTick',[round(freq(4),2) 1 round(freq(3),2)],'XColor','none','XLim',[xunten xoben])
ylim('padded')

nexttile([4 1]);
scatter(PlotRe2,PlotIm2,1,'k')
set(gca,'XTick',[],'YTick',[round(freq(6),2) 0 round(freq(5),2)],'XColor','none','XLim',[xunten xoben])
ylim('padded')

nexttile([4 1]);
scatter(PlotRe3,PlotIm3,1,'k')
set(gca,'XTick',[],'YTick',[round(freq(8),2) -1 round(freq(7),2)],'XColor','none','XLim',[xunten xoben])
ylim('padded')

nexttile([2 1]);
scatter(CharExRe(:,2),CharExIm(:,2),0.8,'k')
set(gca,'XTick',[],'YTick',[-2 round(freq(9),2)],'XColor','none','XLim',[xunten xoben])
ylim('padded')

nexttile([2 1]);
scatter(CharExRe(:,1),CharExIm(:,1),0.8,'k')
set(gca,'YTick',[-3 round(freq(10),2)],'XLim',[xunten xoben])
ylim('padded')

xlabel(tl,'${Re(s_R)}$','interpreter','latex','FontSize', 14);
ylabel(tl,'${Im(s_R)}$','interpreter','latex','FontSize', 14);


saveas(gcf,[pwd '/Plots/Eigenwerte.fig'])
saveas(gcf,[pwd '/Plots/Eigenwerte.png'])


%% 3-Blatt-Rotor Lagerlos gamma1 = 15,238

load("Workspace.mat")
close all


figure('Position', [100 100 700 700]);
s = scatter(damp,freq,100,'linewidth',2);
m = s.Marker;
s.Marker = 'x';
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
set(gca,'XLim',[-2*nu0 2.5*nu0],'YLim',[-2.1*nu0 2.1*nu0])
viscircles([0 0],nu0,'color','k','linewidth',0.8)
line([0 0.75],[0 sqrt(nu0^2-0.75^2)])
xlabel('${Re(s_R)}$','interpreter','latex','FontSize', 14);
ylabel('${Im(s_R)}$','interpreter','latex','FontSize', 14);
daspect([1 1 1])
saveas(gcf,[pwd '/Plots/Schwebeflug.fig'])
saveas(gcf,[pwd '/Plots/Schwebeflug.png'])

%Diagrammgrenzen
xunten = min(CharExRe(length(CharExRe),AnzGl/2))-0.2;
xoben = max(CharExRe(length(CharExRe),AnzGl))+0.2;

xachse = MuMin:SW:MuMax;


% Plotten des Realteils
figure;
plot(xachse,CharExRe,'LineWidth',1.5,'Color','k');
set(gca,'ydir', 'reverse');
xlabel('$\mu$','interpreter','latex','FontSize', 14);
ylabel('${Re(s_R)}$','interpreter','latex','FontSize', 14);
saveas(gcf,[pwd '/Plots/Realteile.fig'])
saveas(gcf,[pwd '/Plots/Realteile.png'])


% Plotten des Imaginaerteils
figure('Position', [10 10 600 300]);
plot(xachse,CharExIm(:,2),'LineWidth',1.5,'Color','k');
ylim('padded')
xlabel('$\mu$','interpreter','latex','FontSize', 14);
ylabel('${Im(s_R)}$','interpreter','latex','FontSize', 14);
saveas(gcf,[pwd '/Plots/Imaginaerteile.fig'])
saveas(gcf,[pwd '/Plots/Imaginaerteile.png'])


PlotRe1 = [CharExRe(:,4:5),CharExRe(:,2:3)];
PlotIm1 = [CharExIm(:,4:5),CharExIm(:,2:3)];

figure('Position', [10 10 600 AnzGl*2*50]);
tl = tiledlayout(AnzGl*2,1);

nexttile([2 1]);
scatter(CharExRe(:,6),CharExIm(:,6),1,'k')
set(gca,'XTick',[],'YTick',[1.5 round(freq(1),2) round(max(CharExIm(:,6)),2)],'XColor','none','XLim',[xunten xoben])
ylim('padded')

nexttile([8 1]);
scatter(PlotRe1,PlotIm1,1,'k')
set(gca,'XTick',[],'YTick',[-1 round(freq(5),2) -0.5 round(freq(4),2) 0 round(freq(3),2) 0.5 round(freq(2),2) 1],'XColor','none','XLim',[xunten xoben])
ylim('padded')


nexttile([2 1]);
scatter(CharExRe(:,1),CharExIm(:,1),1,'k')
set(gca,'XLim',[xunten xoben],'YTick',[round(min(CharExIm(:,1)),2) round(freq(6),2) -1.5])
ylim('padded')

xlabel(tl,'${Re(s_R)}$','interpreter','latex','FontSize', 14);
ylabel(tl,'${Im(s_R)}$','interpreter','latex','FontSize', 14);


saveas(gcf,[pwd '/Plots/Eigenwerte.fig'])
saveas(gcf,[pwd '/Plots/Eigenwerte.png'])

%% 3-Blatt-Rotor Lagerlos gamma2 = 4

load("Workspace.mat")
close all

figure('Position', [100 100 700 700]);
s = scatter(damp,freq,100,'linewidth',2);
m = s.Marker;
s.Marker = 'x';
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
set(gca,'XLim',[-2*nu0 2.5*nu0],'YLim',[-2.1*nu0 2.1*nu0])
viscircles([0 0],nu0,'color','k','linewidth',0.8)
line([0 0.75],[0 sqrt(nu0^2-0.75^2)])
xlabel('${Re(s_R)}$','interpreter','latex','FontSize', 14);
ylabel('${Im(s_R)}$','interpreter','latex','FontSize', 14);
daspect([1 1 1])
saveas(gcf,[pwd '/Plots/Schwebeflug.fig'])
saveas(gcf,[pwd '/Plots/Schwebeflug.png'])

%Diagrammgrenzen
xunten = min(CharExRe(length(CharExRe),AnzGl))-0.2;
xoben = max(CharExRe(length(CharExRe),1))+0.2;

xachse = MuMin:SW:MuMax;


% Plotten des Realteils
figure;
plot(xachse,CharExRe,'LineWidth',1.5,'Color','k');
set(gca,'ydir', 'reverse','YLim',[xunten xoben]);
xlabel('$\mu$','interpreter','latex','FontSize', 14);
ylabel('${Re(s_R)}$','interpreter','latex','FontSize', 14);
saveas(gcf,[pwd '/Plots/Realteile.fig'])
saveas(gcf,[pwd '/Plots/Realteile.png'])


% Plotten des Imaginaerteils
figure('Position', [10 10 600 300]);
plot(xachse,CharExIm(:,2),'LineWidth',1.5,'Color','k');
ylim('padded')
xlabel('$\mu$','interpreter','latex','FontSize', 14);
ylabel('${Im(s_R)}$','interpreter','latex','FontSize', 14);
saveas(gcf,[pwd '/Plots/Imaginaerteile.fig'])
saveas(gcf,[pwd '/Plots/Imaginaerteile.png'])


PlotRe1 = [CharExRe(:,4),CharExRe(:,3)];
PlotIm1 = [CharExIm(:,4),CharExIm(:,3)];


figure('Position', [10 10 600 AnzGl*2*50]);
tl = tiledlayout(AnzGl*2,1);

nexttile([2 1]);
scatter(CharExRe(:,1),CharExIm(:,1),1,'k')
set(gca,'XTick',[],'YTick',[2 round(freq(1),2)],'XColor','none','XLim',[xunten xoben])
ylim('padded')

nexttile([2 1]);
scatter(CharExRe(:,2),CharExIm(:,2),1,'k')
set(gca,'XTick',[],'YTick',[1 round(freq(2),2)],'XColor','none','XLim',[xunten xoben])
ylim('padded')

nexttile([4 1]);
scatter(PlotRe1,PlotIm1,1,'k')
set(gca,'XTick',[],'YTick',[round(freq(4),2) 0 round(freq(3),2)],'XColor','none','XLim',[xunten xoben])
ylim('padded')

nexttile([2 1]);
scatter(CharExRe(:,5),CharExIm(:,5),1,'k')
set(gca,'XLim',[xunten xoben],'YTick',[round(freq(5),2) -1],'XColor','none')
ylim('padded')

nexttile([2 1]);
scatter(CharExRe(:,6),CharExIm(:,6),1,'k')
set(gca,'YTick',[round(freq(6),2) -2],'XLim',[xunten xoben])
ylim('padded')

xlabel(tl,'${Re(s_R)}$','interpreter','latex','FontSize', 14);
ylabel(tl,'${Im(s_R)}$','interpreter','latex','FontSize', 14);


saveas(gcf,[pwd '/Plots/Eigenwerte.fig'])
saveas(gcf,[pwd '/Plots/Eigenwerte.png'])


%% 4-Blatt-Rotor Lagerlos gamma1 = 8,16

load("Workspace.mat")
close all


%Diagrammgrenzen
xunten = min(CharExRe(:,AnzGl))-0.2;
xoben = max(CharExRe(:,1))+0.2;

figure('Position', [100 100 700 700]);
s = scatter(damp,freq,100,'linewidth',2);
m = s.Marker;
s.Marker = 'x';
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
set(gca,'XLim',[-2*nu0 2.5*nu0],'YLim',[-2.1*nu0 2.1*nu0])
viscircles([0 0],nu0,'color','k','linewidth',0.8)
line([0 0.75],[0 sqrt(nu0^2-0.75^2)])
xlabel('${Re(s_R)}$','interpreter','latex','FontSize', 14);
ylabel('${Im(s_R)}$','interpreter','latex','FontSize', 14);
daspect([1 1 1])
saveas(gcf,[pwd '/Plots/Schwebeflug.fig'])
saveas(gcf,[pwd '/Plots/Schwebeflug.png'])

xachse = MuMin:SW:MuMax;


% Plotten des Realteils
figure;
plot(xachse,CharExRe,'LineWidth',1.5,'Color','k');
set(gca,'ydir', 'reverse','YLim',[xunten xoben]);
xlabel('$\mu$','interpreter','latex','FontSize', 14);
ylabel('${Re(s_R)}$','interpreter','latex','FontSize', 14);
saveas(gcf,[pwd '/Plots/Realteile.fig'])
saveas(gcf,[pwd '/Plots/Realteile.png'])


% Plotten des Imaginaerteils
figure('Position', [10 10 600 300]);
plot(xachse,CharExIm(:,2),'LineWidth',1.5,'Color','k');
%scatter(xachse,CharExIm,1)
ylim('padded')
xlabel('$\mu$','interpreter','latex','FontSize', 14);
ylabel('${Im(s_R)}$','interpreter','latex','FontSize', 14);
saveas(gcf,[pwd '/Plots/Imaginaerteile.fig'])
saveas(gcf,[pwd '/Plots/Imaginaerteile.png'])


PlotRe1 = [CharExRe(:,4),CharExRe(:,5)];
PlotIm1 = [CharExIm(:,4),CharExIm(:,5)];


figure('Position', [10 10 600 6*2*50]);
tl = tiledlayout(6*2,1);

nexttile([2 1]);
scatter(CharExRe(:,1),CharExIm(:,1),0.8,'k')
set(gca,'XTick',[],'YTick',[2 round(max(CharExIm(:,1)),2)],'XColor','none','XLim',[xunten xoben])
ylim('padded')

nexttile([2 1]);
scatter(CharExRe(:,2),CharExIm(:,2),0.8,'k')
set(gca,'XTick',[],'YTick',[1 round(max(CharExIm(:,2)),2)],'XColor','none','XLim',[xunten xoben])
ylim('padded')

nexttile([4 1]);
scatter(PlotRe1,PlotIm1,1,'k')
set(gca,'XTick',[],'YTick',[round(min(PlotIm1(:,2)),2) 0 round(max(PlotIm1(:,1)),2)],'XColor','none','XLim',[xunten xoben])
ylim('padded')

nexttile([2 1]);
scatter(CharExRe(:,7),CharExIm(:,7),0.8,'k')
set(gca,'XTick',[],'YTick',[round(min(CharExIm(:,7)),2) -1],'XColor','none','XLim',[xunten xoben])
ylim('padded')

nexttile([2 1]);
scatter(CharExRe(:,8),CharExIm(:,8),0.8,'k')
set(gca,'XLim',[xunten xoben],'YTick',[round(min(CharExIm(:,8)),2) -1])
ylim('padded')

xlabel(tl,'${Re(s_R)}$','interpreter','latex','FontSize', 14);
ylabel(tl,'${Im(s_R)}$','interpreter','latex','FontSize', 14);


saveas(gcf,[pwd '/Plots/Eigenwerte.fig'])
saveas(gcf,[pwd '/Plots/Eigenwerte.png'])


%% 4-Blatt-Rotor Lagerlos gamma2 = 4

load("Workspace.mat")
close all

figure('Position', [100 100 700 700]);
s = scatter(damp,freq,100,'linewidth',2);
m = s.Marker;
s.Marker = 'x';
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
set(gca,'XLim',[-2*nu0 2.5*nu0],'YLim',[-2.1*nu0 2.1*nu0])
viscircles([0 0],nu0,'color','k','linewidth',0.8)
line([0 0.75],[0 sqrt(nu0^2-0.75^2)])
xlabel('${Re(s_R)}$','interpreter','latex','FontSize', 14);
ylabel('${Im(s_R)}$','interpreter','latex','FontSize', 14);
daspect([1 1 1])
saveas(gcf,[pwd '/Plots/Schwebeflug.fig'])
saveas(gcf,[pwd '/Plots/Schwebeflug.png'])

%Diagrammgrenzen
xunten = min(CharExRe(:,AnzGl))-0.2;
xoben = max(CharExRe(:,1))+0.2;

xachse = MuMin:SW:MuMax;

% Plotten des Realteils
figure;
plot(xachse,CharExRe,'LineWidth',1.5,'Color','k');
set(gca,'ydir', 'reverse','YLim',[xunten xoben]);
xlabel('$\mu$','interpreter','latex','FontSize', 14);
ylabel('${Re(s_R)}$','interpreter','latex','FontSize', 14);
saveas(gcf,[pwd '/Plots/Realteile.fig'])
saveas(gcf,[pwd '/Plots/Realteile.png'])


% Plotten des Imaginaerteils
figure('Position', [10 10 600 300]);
plot(xachse,CharExIm(:,2),'LineWidth',1.5,'Color','k');
%scatter(xachse,CharExIm,1)
ylim('padded')
xlabel('$\mu$','interpreter','latex','FontSize', 14);
ylabel('${Im(s_R)}$','interpreter','latex','FontSize', 14);
saveas(gcf,[pwd '/Plots/Imaginaerteile.fig'])
saveas(gcf,[pwd '/Plots/Imaginaerteile.png'])


PlotRe1 = [CharExRe(:,4),CharExRe(:,5)];
PlotIm1 = [CharExIm(:,4),CharExIm(:,5)];


figure('Position', [10 10 600 6*2*50]);
tl = tiledlayout(6*2,1);

nexttile([2 1]);
scatter(CharExRe(:,1),CharExIm(:,1),0.8,'k')
set(gca,'XTick',[],'YTick',[2 round(max(CharExIm(:,1)),2)],'XColor','none','XLim',[xunten xoben])
ylim('padded')

nexttile([2 1]);
scatter(CharExRe(:,2),CharExIm(:,2),0.8,'k')
set(gca,'XTick',[],'YTick',[1 round(max(CharExIm(:,2)),2)],'XColor','none','XLim',[xunten xoben])
ylim('padded')

nexttile([4 1]);
scatter(PlotRe1,PlotIm1,1,'k')
set(gca,'XTick',[],'YTick',[round(min(PlotIm1(:,2)),2) 0 round(max(PlotIm1(:,1)),2)],'XColor','none','XLim',[xunten xoben])
ylim('padded')

nexttile([2 1]);
scatter(CharExRe(:,7),CharExIm(:,7),0.8,'k')
set(gca,'XTick',[],'YTick',[round(min(CharExIm(:,7)),2) -1],'XColor','none','XLim',[xunten xoben])
ylim('padded')

nexttile([2 1]);
scatter(CharExRe(:,8),CharExIm(:,8),0.8,'k')
set(gca,'XLim',[xunten xoben],'YTick',[round(min(CharExIm(:,8)),2) -1])
ylim('padded')

xlabel(tl,'${Re(s_R)}$','interpreter','latex','FontSize', 14);
ylabel(tl,'${Im(s_R)}$','interpreter','latex','FontSize', 14);


saveas(gcf,[pwd '/Plots/Eigenwerte.fig'])
saveas(gcf,[pwd '/Plots/Eigenwerte.png'])



%% 4-Blatt-Lagerlos gamma1 = 8,16, konstante Koeffizienten

load("Workspace.mat")
close all

%Diagrammgrenzen
xunten = min(CharExRe(:,AnzGl))-0.2;
xoben = max(CharExRe(:,1))+0.2;

figure('Position', [100 100 700 700]);
s = scatter(damp,freq,100,'linewidth',2);
m = s.Marker;
s.Marker = 'x';
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
set(gca,'XLim',[-2*nu0 2.5*nu0],'YLim',[-2.1*nu0 2.1*nu0])
viscircles([0 0],nu0,'color','k','linewidth',0.8)
line([0 0.75],[0 sqrt(nu0^2-0.75^2)])
xlabel('${Re(s_R)}$','interpreter','latex','FontSize', 14);
ylabel('${Im(s_R)}$','interpreter','latex','FontSize', 14);
daspect([1 1 1])
saveas(gcf,[pwd '/Plots/Schwebeflug.fig'])
saveas(gcf,[pwd '/Plots/Schwebeflug.png'])

xachse = MuMin:SW:MuMax;


% Plotten des Realteils
figure;
plot(xachse,CharExRe,'LineWidth',1.5,'Color','k');
set(gca,'ydir', 'reverse','YLim',[xunten xoben]);
xlabel('$\mu$','interpreter','latex','FontSize', 14);
ylabel('${Re(s_R)}$','interpreter','latex','FontSize', 14);
saveas(gcf,[pwd '/Plots/Realteile.fig'])
saveas(gcf,[pwd '/Plots/Realteile.png'])


% Plotten des Imaginaerteils
figure('Position', [10 10 600 300]);
plot(xachse,CharExIm(:,2),'LineWidth',1.5,'Color','k');
%scatter(xachse,CharExIm,1)
ylim('padded')
xlabel('$\mu$','interpreter','latex','FontSize', 14);
ylabel('${Im(s_R)}$','interpreter','latex','FontSize', 14);
saveas(gcf,[pwd '/Plots/Imaginaerteile.fig'])
saveas(gcf,[pwd '/Plots/Imaginaerteile.png'])


PlotRe1 = [CharExRe(:,4),CharExRe(:,5)];
PlotIm1 = [CharExIm(:,4),CharExIm(:,5)];


figure('Position', [10 10 600 6*2*50]);
tl = tiledlayout(6*2,1);

nexttile([2 1]);
scatter(CharExRe(:,1),CharExIm(:,1),0.8,'k')
set(gca,'XTick',[],'YTick',[2 round(max(CharExIm(:,1)),2)],'XColor','none','XLim',[xunten xoben])
ylim('padded')

nexttile([2 1]);
scatter(CharExRe(:,2),CharExIm(:,2),0.8,'k')
set(gca,'XTick',[],'YTick',[1 round(max(CharExIm(:,2)),2)],'XColor','none','XLim',[xunten xoben])
ylim('padded')

nexttile([4 1]);
scatter(PlotRe1,PlotIm1,1,'k')
set(gca,'XTick',[],'YTick',[round(min(PlotIm1(:,2)),2) 0 round(max(PlotIm1(:,1)),2)],'XColor','none','XLim',[xunten xoben])
ylim('padded')

nexttile([2 1]);
scatter(CharExRe(:,7),CharExIm(:,7),0.8,'k')
set(gca,'XTick',[],'YTick',[round(min(CharExIm(:,7)),2) -1],'XColor','none','XLim',[xunten xoben])
ylim('padded')

nexttile([2 1]);
scatter(CharExRe(:,8),CharExIm(:,8),0.8,'k')
set(gca,'XLim',[xunten xoben],'YTick',[round(min(CharExIm(:,8)),2) -1])
ylim('padded')


xlabel(tl,'${Re(s_R)}$','interpreter','latex','FontSize', 14);
ylabel(tl,'${Im(s_R)}$','interpreter','latex','FontSize', 14);

saveas(gcf,[pwd '/Plots/Eigenwerte.fig'])
saveas(gcf,[pwd '/Plots/Eigenwerte.png'])

%% Einzelblattkoordinaten

load("Workspace.mat")
close all

%Diagrammgrenzen
xunten = min(CharExRe(length(CharExRe),2))-0.2;
xoben = max(CharExRe(length(CharExRe),1))+0.2;

xachse = MuMin:SW:MuMax;

figure('Position', [100 100 700 700]);
s = scatter(damp,freq,100,'linewidth',2);
m = s.Marker;
s.Marker = 'x';
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
set(gca,'XLim',[-2*nu0 2.5*nu0],'YLim',[-2.1*nu0 2.1*nu0])
viscircles([0 0],nu0,'color','k','linewidth',0.8)
line([0 0.75],[0 sqrt(nu0^2-0.75^2)])
xlabel('${Re(s_R)}$','interpreter','latex','FontSize', 14);
ylabel('${Im(s_R)}$','interpreter','latex','FontSize', 14);
daspect([1 1 1])
saveas(gcf,[pwd '/Plots/Schwebeflug.fig'])
saveas(gcf,[pwd '/Plots/Schwebeflug.png'])

% Plotten des Realteils
figure;
plot(xachse,CharExRe,'LineWidth',1.5,'Color','k');
set(gca,'ydir', 'reverse');
xlabel('$\mu$','interpreter','latex','FontSize', 14);
ylabel('${Re(s_R)}$','interpreter','latex','FontSize', 14);
saveas(gcf,[pwd '/Plots/Realteile.fig'])
saveas(gcf,[pwd '/Plots/Realteile.png'])

% Plotten des Imaginaerteils
figure('Position', [10 10 600 300]);
plot(xachse,CharExIm(:,1),'LineWidth',1.5,'Color','k');
ylim('padded')
xlabel('$\mu$','interpreter','latex','FontSize', 14);
ylabel('${Im(s_R)}$','interpreter','latex','FontSize', 14);
saveas(gcf,[pwd '/Plots/Imaginaerteile.fig'])
saveas(gcf,[pwd '/Plots/Imaginaerteile.png'])

% Plotten der Eigenwerte in der komplexen Ebene
figure('Position', [10 10 600 AnzGl*2*50]);
tl = tiledlayout(AnzGl*2,1);

nexttile([2 1]);
scatter(CharExRe(:,1),CharExIm(:,1),1,'k')
set(gca,'XTick',[],'YTick',[1 round(freq(1),2)],'XColor','none','XLim',[xunten xoben])
ylim('padded')

nexttile([2 1]);
scatter(CharExRe(:,2),CharExIm(:,2),1,'k')
set(gca,'YTick',[round(freq(2),2) -1],'XLim',[xunten xoben])
ylim('padded')

xlabel(tl,'${Re(s_R)}$','interpreter','latex','FontSize', 14);
ylabel(tl,'${Im(s_R)}$','interpreter','latex','FontSize', 14);

saveas(gcf,[pwd '/Plots/Eigenwerte.fig'])
saveas(gcf,[pwd '/Plots/Eigenwerte.png'])

%% Ende
close all