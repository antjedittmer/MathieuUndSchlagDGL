%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Erzeugen der Diagramme fuer die Verlaeufe der Transitionsmatrizen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MathieuDiagrammeTransitionsmatrizen_GeschwX1aufYAchse
clc; close all;

isInterpTex = 0; %1: Tex interpreter, 0: Latex interpreter
plot2D = 0; %1: 2D und 3D plots, 0: nur 3D plots erstellen
useAspect = 0;

if isInterpTex == 1 || verLessThan('matlab', '9.8')
    isInterpTex = 1; % Bei R2007 wird der Tex-Interpreter eingeschaltet
    strInterp = 'tex';
else
    isInterpTex = 0;
    strInterp = 'latex';
end

nu_02 = 5;
nu_C2 = 5;

t0 = 0.0;
T = 2*pi;
tspan = t0:0.0001:T;
Diagonal = diag(ones(2,1));


fDir = 'figureFolder_GeschwX1aufYAchse'; % Ordner Abbildungen
if ~isdir(fDir) %#ok<ISDIR>
    mkdir(fDir)
end

fDirAspect = 'figureFolderAspectPhi';
if ~isdir(fDirAspect) %#ok<ISDIR>
    mkdir(fDirAspect)
end

% D = [0, 0.1, 0.15, 0.8, 1];
% D = 0.1;
DVec = [0.15; 0.001; 0.2];
strDVec = {sprintf('%2.2f',DVec(1)); sprintf('%2.3f',DVec(2)); ...
    sprintf('%2.1f',DVec(3))};

% Für die grafische Darstellung
strP.hi = '$\phi \; [-]$';
strP.hiDot = '$\stackrel{\ast}{\phi} \; [-]$';
strP.si = '$\psi \; [rad]$';
strP.ystr = 'Zustandsgrößen';

if isInterpTex == 1
    strP.hi = '\phi [-]';
    strP.hiDot = '\phi^* [-]';
    strP.si = '\psi [rad]';
end


noF = 1;
for dIdx = 1: length(DVec)
    D = DVec(dIdx);
    strD = strDVec{dIdx};

    sol1 = ode45(@(psi,x)MathieuDGL(psi,x,D,nu_02,nu_C2),[t0,T],Diagonal(:,1));
    sol2 = ode45(@(psi,x)MathieuDGL(psi,x,D,nu_02,nu_C2),[t0,T],Diagonal(:,2));

    y1 = deval(sol1,tspan);
    y2 = deval(sol2,tspan);

    Pos1 = y1(2,:)' ;
    Gesch1 = y1(1,:)' ;

    Pos2 = y2(2,:)' ;
    Gesch2 = y2(1,:)';

    %% Graphiken: 2D und 3D, 1. und 2. Spaltenvektor: 1. Spaltenvektor
    AxisLimitsYZ.Y1 = [-2 2]; %[-1,1];
    AxisLimitsYZ.Z1 = [-1,1]; %[-2 2];
    labelPos.x = [7.6 0.2 0];
    labelPos.y = [0 1.2 0.1];
    labelPos.z = [0.7 0 0.7]; %[0.7 0 1.3];

    titleStr = strrep(sprintf('%s,%s','1. Spaltenvektor der Monodromiematrix fuer e = [1;0]; D = ',strD),',','');
    pngname = strrep(sprintf('MathieuDiagram_D%2.1e_x1_2D',D),'.','dot');

    if plot2D == 1
        pngfile = fullfile(fDir,pngname);
        noF = plot3D(noF,tspan,Gesch1,Pos1,T,AxisLimitsYZ,strP,labelPos,titleStr,pngfile,plot2D,isInterpTex);
    end

    pngname3D = strrep(pngname,'2D','3D');
    pngfile = fullfile(fDir,pngname3D);
    [noF,hAxis] = plot3D(noF,tspan,Gesch1,Pos1,T,AxisLimitsYZ,strP,labelPos, titleStr,pngfile,0,isInterpTex);

    if useAspect == 1
        daspect([1 1 1]);
        if isInterpTex == 0
            ylabel(hAxis,strP.hi,'interpreter',strInterp,'Position',labelPos.y.*[1,1.2,1])
        end
        view([60,15])
        pngfileAspect = strrep(pngfile,fDir,fDirAspect);
        print(pngfileAspect, '-dpng')
    end

    %% Graphiken: 2D und 3D, 1. und 2. Spaltenvektor: 1. Spaltenvektor
    AxisLimitsYZ.Y1 = [-4.5 5.1]; %[-2,2.4];
    AxisLimitsYZ.Z1 = [-2,2.4]; %[-4.5 5.1];
    multY = 2;
    if D < 0.1
        labelPos.x = [8.5 0.2 0];
        multY = 1.5;
    end
    labelPos.y = [0 2.4 0.8];
    labelPos.z = [0 1 2]; %[0 0.5 2.8];

    titleStr2 = strrep(strrep(titleStr,'1.','2.'),'[1;0]', '[0;1]');
    pngname2 = strrep(pngname,'x1','x2');

    if plot2D == 1
        pngfile = fullfile(fDir,pngname2);
        noF = plot3D(noF,tspan,Gesch2,Pos2,T,AxisLimitsYZ,strP,labelPos,titleStr2,pngfile,plot2D,isInterpTex);
    end

    pngname2_3D = strrep(pngname2,'2D','3D');
    pngfile = fullfile(fDir,pngname2_3D);
    [noF,hAxis] = plot3D(noF,tspan,Gesch2,Pos2,T,AxisLimitsYZ,strP,labelPos,titleStr2,pngfile,0,isInterpTex);

    if useAspect == 1
        daspect([1 1 1]);
        if isInterpTex == 0
            labelPos.x = [8.5 0.2 0];
            xlabel(hAxis,strP.si,'interpreter',strInterp,'Position',labelPos.x)
            ylabel(hAxis,strP.hi,'interpreter',strInterp,'Position',labelPos.y .* [1,multY,0.1])
            zlabel(hAxis,strP.hiDot,'interpreter',strInterp,'Position',labelPos.z .* [1,2,1],'Rotation',0)
        end
        view([40,15])
        pngfileAspect = strrep(pngfile,fDir,fDirAspect);
        print(pngfileAspect, '-dpng')
    end

end

end

function [noF,hAxis]  = plot3D(noF,tspan,Gesch1,Pos1,T,AxisLimitsYZ,strP,labelPos,titleStr,pngfile,plot2D,isInterpTex)

if nargin < 11 || isempty(plot2D)
    plot2D = 0;
end

strInterp = 'tex';
if nargin <12 || isInterpTex == 0
    isInterpTex = 0;
    strInterp = 'latex';
end

figure(noF); noF = noF + 1; % 3D, 1. Spaltenvektor

if plot2D % 2D: Default x,y-label, Legende
    plot(tspan,Gesch1,tspan,Pos1,'--')
    xlabel(strP.si,'interpreter',strInterp);
    ylabel(strP.ystr);
    axis tight; grid on;
    % legend(strP.hi,strP.hiDot,'interpreter',strInterp,'location','southeast')
    legend(strP.hiDot,strP.hi,'interpreter',strInterp,'location','southeast')
else
    plot3(tspan,Gesch1,Pos1)
    lim2 = max(max(abs(Gesch1)),1);
    axis([0 T -lim2 lim2])
end
title(titleStr); % title

if isInterpTex
    hAxis = get(gca);
else
    hAxis = gca;
end
hAxis.XTick = 0:pi/2:2*pi;
hAxis.XTickLabel = {'0','\pi/2','\pi','3\pi/2','2\pi'};

if plot2D
    if isInterpTex == 1
        hAxis.XAxisLocation = 'origin';
        hAxis.YAxisLocation = 'origin';
    end
else
    if isInterpTex == 0 %only use this if the interpreter is Latex
        hAxis.XRuler.FirstCrossoverValue  = 0; % X crossover with Y axis
        hAxis.YRuler.FirstCrossoverValue  = 0; % Y crossover with X axis
        hAxis.ZRuler.FirstCrossoverValue  = 0; % Z crossover with X axis
        hAxis.ZRuler.SecondCrossoverValue = 0; % Z crossover with Y axis
        hAxis.XRuler.SecondCrossoverValue = 0; % X crossover with Z axis
        hAxis.YRuler.SecondCrossoverValue = 0; % Y crossover with Z axis
        xlabel(hAxis,strP.si,'interpreter',strInterp,'Position',labelPos.x)
        ylabel(hAxis,strP.hiDot,'interpreter',strInterp,'Position',labelPos.y)
        zlabel(hAxis,strP.hi,'interpreter',strInterp,'Position',labelPos.z,'Rotation',0)
    else
        xlabel(strP.si)
        ylabel(strP.hiDot)
        zlabel(strP.hiDot)
        grid on;
    end
    hAxis.YAxis.Limits = AxisLimitsYZ.Y1;
    hAxis.ZAxis.Limits = AxisLimitsYZ.Z1;

% view([130 10])
  view([40 20])
end

print(pngfile, '-dpng')
end