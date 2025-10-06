%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation of the Strutt Charts within the limits of nu_02 and nu_C2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc ; clear all ; close all;
% Step size
SW = 0.1;
% Limits for nu_02 and nu_C2
unt0 = 0; % lower limit for nu_02
untC = 0; % lower limit for nu_C2
ob0 = 9;  % upper limit for nu_02
obC = 9;  % upper limit for nu_C2
% Counter variables
l=1;
o=1;
% Number of equations in the DGL system (Differential Equation System)
Nz=2;
% Parameters
D = 0.001;
t0 = 0.0;
T = 2*pi;
tspan = t0:0.0001:T;
% Initial conditions
Diagonal=diag(ones(Nz,1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preparation of arrays to be filled in the loops for program acceleration
Monodromy = zeros(Nz);
CharMult =  zeros(length(unt0:SW:ob0),2);
CharEx = zeros(length(unt0:SW:ob0),Nz*2);
plotwertstabil = zeros(length(unt0:SW:ob0)*length(untC:SW:obC),3);
plotwertinstabil = zeros(length(unt0:SW:ob0)*length(untC:SW:obC),3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for m = unt0:SW:ob0
    for n = untC:SW:obC   
        nu_02 = m ;
        nu_C2 = n ;
        
        % The loop numerically integrates Nz times consecutively for the
        % column vectors of the identity matrix. The solution evaluated at T
        % is written into MonoVek, and these vectors are assembled into the
        % Monodromy Matrix.
        for k=1:Nz
            sol = ode45(@(psi,x)MathieuDGL(psi,x,D,nu_02,nu_C2),[t0,T],Diagonal(:,k));
            MonoVek = deval(sol,T);
            Monodromie(:,k) = MonoVek;
        end
    
        % characteristic multipliers (eigenvalues of the Monodromy Matrix)
        p = [1 -(Monodromie(1,1)+Monodromie(2,2)) (Monodromie(1,1)*Monodromie(2,2))-Monodromie(1,2)*Monodromie(2,1)];
        e = roots(p).';
        CharMult(l,:) = e;
        % Absolute value of the eigenvalues for stability analysis
        eBetr = abs(e);
        
        Real = 1/T * log(e);
        
        % characteristic exponents
        if nu_C2==nu_02
            CharEx(o,:) = [nu_02,nu_C2,Real(1,1),Real(1,2)];
            o=o+1;
        end
        
        % 1=stable , -1=unstable Note do not forget (to =1)
        % sorting the parameter combinations into an array with stable and
        % with unstable combinations of nu_02 and nu_C2 based on the
        % characteristic multiplier
        b=0;
        if eBetr(1)<1 && eBetr(2)<1
            b=1;
            plotwertstabil(l,:) = [nu_02,nu_C2,b];
        elseif eBetr(1)>1 || eBetr(2)>1
            b=-1;
            plotwertinstabil(l,:) = [nu_02,nu_C2,b];
        end
    
        l=l+1;  
    end
end
% Deleting the purely 0-rows
plotwertstabil = plotwertstabil(any(plotwertstabil,2),:);
plotwertinstabil = plotwertinstabil(any(plotwertinstabil,2),:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Graphical Representation: Strutt Chart (Plotting all pairs of values
% nu_02/nu_C2 that lead to a stable solution) as well as Real parts of the
% characteristic exponents
figure('Position', [10 10 600 600])
tiledlayout(9,1);
nexttile([7 1])
% Strutt Chart
scatter(plotwertstabil(:,1),plotwertstabil(:,2),1)
ylabel('$\rm{Parameter} \; \nu_C^2 \;\;\rm{[-]}$','interpreter','latex','Position', [-0.5 nu_C2/2],'FontSize', 10);
%xlabel('$\rm{Parameter} \; \nu_0^2 \;\;\rm{[-]}$','interpreter','latex','FontSize', 10);
% Real parts char. Exponents
nexttile([2 1])
xachse = unt0:SW:ob0;
plot(xachse,CharEx(:,3:4))
xlabel('$\rm{Parameter} \; \nu_0^2 \;\;\rm{[-]}$','interpreter','latex','FontSize', 10);
ylabel('$Re(s_R) \;\; \rm{[-]}$','interpreter','latex','Position', [-0.5 -D],'FontSize', 10);