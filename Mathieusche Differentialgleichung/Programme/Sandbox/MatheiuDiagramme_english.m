%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generation of the diagrams for the time histories of the transition matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;
nu_02 = 5;
nu_C2 = 5;
t0 = 0.0;
T = 2*pi;
tspan = t0:0.0001:T;
Diagonal = diag(ones(2,1));
plot2D = 0;
fDir = 'figureFolder4'; % Folder for figures (Abbildungen)
if ~isdir(fDir) %#ok<ISDIR>
    mkdir(fDir)
end
% D = [0, 0.1, 0.15, 0.8, 1];
% D = 0.1;
DVec = [0.15; 0.001; 0.2];
strDVec = {sprintf('%2.2f',DVec(1)); sprintf('%2.3f',DVec(2)); ...
    sprintf('%2.1f',DVec(3))};
% strDVec = cell(length(DVec),1);
% for idx = 1: length(DVec)
%     strDVec{idx} = sprintf('%2.2f',DVec(idx));
% end
noF = 1;
for dIdx = 1: length(DVec)
    D = DVec(dIdx);
    strD = strDVec{dIdx};
    % graphical representation (grafische Darstellung)
    sol1 = ode45(@(psi,x)MathieuDGL(psi,x,D,nu_02,nu_C2),[t0,T],Diagonal(:,1));
    sol2 = ode45(@(psi,x)MathieuDGL(psi,x,D,nu_02,nu_C2),[t0,T],Diagonal(:,2));
    y1 = deval(sol1,tspan);
    y2 = deval(sol2,tspan);
    Pos1 = y1(1,:)' ;
    Gesch1 = y1(2,:)' ;
    Pos2 = y2(1,:)' ;
    Gesch2 = y2(2,:)';
    strPhi = '$\phi \; [-]$';
    strPhiDot = '$\stackrel{\ast}{\phi} \; [1/s]$';
    strPsi = '$\psi \; [rad]$';
    %% Graphics: 2D and 3D, 1st and 2nd column vectors
    ystr = 'State variables (Zustandsgrößen)';
    titleStr = strrep(sprintf('%s,%s','1st column vector of the Monodromy Matrix for e = [1;0]; D = ',strD),',','');
    pngname = strrep(sprintf('MathieuDiagram_D%2.1e_x1_2D',D),'.','dot');
    if plot2D == 1
        figure(noF); noF = noF + 1; % 2D, 1st column vector
        plot(tspan,Pos1,tspan,Gesch1,'--')
        xlabel('$\psi$','interpreter','latex');
        ylabel(ystr);
        axis tight;
        ax = gca;
        ax.XAxisLocation = 'origin';
        ax.YAxisLocation = 'origin';
        set(gca,'XTick',0:pi/2:2*pi)
        set(gca,'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'})
        legend(strPhi,strPhiDot,'interpreter','latex','location','southeast')
        grid on;
        title(titleStr)
        pngfile = fullfile(fDir,pngname);
        print(pngfile, '-dpng')
    end
    figure(noF); noF = noF + 1; % 3D, 1st column vector
    plot3(tspan,Gesch1,Pos1)
    axis([0 T -1 1])
    set(gca,'XTick',0:pi/2:2*pi)
    set(gca,'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'})
    hAxis = gca;
    hAxis.XRuler.FirstCrossoverValue  = 0; % X crossover with Y axis
    hAxis.YRuler.FirstCrossoverValue  = 0; % Y crossover with X axis
    hAxis.ZRuler.FirstCrossoverValue  = 0; % Z crossover with X axis
    hAxis.ZRuler.SecondCrossoverValue = 0; % Z crossover with Y axis
    hAxis.XRuler.SecondCrossoverValue = 0; % X crossover with Z axis
    hAxis.YRuler.SecondCrossoverValue = 0; % Y crossover with Z axis
    xlabel(hAxis,strPsi,'interpreter','latex','Position',[7.6 0.2 0])
    ylabel(hAxis,strPhiDot,'interpreter','latex','Position',[0 1.2 0.1])
    zlabel(hAxis,strPhi,'interpreter','latex','Position',[0.7 0 1.3],'Rotation',0)
    if D==0.1 && nu_02 == 5 && nu_C2 == 5
        line( [2*pi 2*pi] , [0 0.189] , [0 0.723] , 'Color','black','LineStyle','--')
        text(6.5,0.2,0.5,'$\underline{\phi}_1(T,0)$','interpreter','latex')
    end
    hAxis.YAxis.Limits = [-1,1];
    hAxis.ZAxis.Limits = [-2 2];
    title(titleStr)
    view([50 20])
    pngname3D = strrep(pngname,'2D','3D');
    pngfile = fullfile(fDir,pngname3D);
    print(pngfile, '-dpng')
    titleStr2 = strrep(strrep(titleStr,'1.','2.'),'[1;0]', '[0;1]');
    pngname2 = strrep(pngname,'x1','x2');
    if plot2D == 1
        figure(noF); noF = noF + 1; % 2D, 1st column vector
        plot(tspan,Pos2,tspan,Gesch2,'--')
        ax = gca;
        ax.XAxisLocation = 'origin';
        ax.YAxisLocation = 'origin';
        set(gca,'XTick',0:pi/2:2*pi)
        set(gca,'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'})
        xlabel(strPsi,'interpreter','latex');
        ylabel(ystr);
        legend(strPhi,strPhiDot,'interpreter','latex','location','southeast')
        grid on;
        title(titleStr2)
        pngfile = fullfile(fDir,pngname2);
        print(pngfile, '-dpng')
    end
    figure(noF); noF = noF + 1;
    plot3(tspan,Gesch2,Pos2)
    axMax = max(max(Gesch2) *1.05, 2);
    aMin = min(min(Gesch2) *1.05, -2);
    axis([0 T aMin axMax])
    set(gca,'XTick',0:pi/2:2*pi)
    set(gca,'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'})
    hAxis = gca;
    hAxis.XRuler.FirstCrossoverValue  = 0; % X crossover with Y axis
    hAxis.YRuler.FirstCrossoverValue  = 0; % Y crossover with X axis
    hAxis.ZRuler.FirstCrossoverValue  = 0; % Z crossover with X axis
    hAxis.ZRuler.SecondCrossoverValue = 0; % Z crossover with Y axis
    hAxis.XRuler.SecondCrossoverValue = 0; % X crossover with Z axis
    hAxis.YRuler.SecondCrossoverValue = 0; % Y crossover with Z axis
    xlabel(hAxis,strPsi,'interpreter','latex','Position',[8.5 0.2 0])
    ylabel(hAxis,strPhiDot,'interpreter','latex','Position',[0 2.4 0.8])
    zlabel(hAxis,strPhi,'interpreter','latex','Position',[0 0.5 2.8],'Rotation',0)
    if D==0.1 && nu_02 == 5 && nu_C2 == 5
        line( [2*pi 2*pi] , [0 0.76] , [0 1.4] , 'Color','black','LineStyle','--')
        text(7,0.2,1,'$\underline{\phi}_2(T,0)$','interpreter','latex')
    end
    titleStr2 = strrep(strrep(titleStr,'1.','2.'),'[1;0]', '[0;1]');
    title(titleStr2)
    view([50 20])
    ayMax = max(max(Pos2) *1.05, 2);
    yMin = min(min(Pos2) *1.05, -2);
    hAxis.YAxis.Limits = [-2,2.4];
    hAxis.ZAxis.Limits = [-4.5 5.1];
    %daspect([4 ayMax ayMax])
    pngname2_3D = strrep(pngname2,'2D','3D');
    pngfile = fullfile(fDir,pngname2_3D);
    print(pngfile, '-dpng')
end