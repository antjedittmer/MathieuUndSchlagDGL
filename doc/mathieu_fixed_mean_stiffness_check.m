% Control experiment for Floquet_ModalTracking_Findings.tex:
% frequency of the dominant harmonic of the Floquet mode of the damped
% Mathieu ODE  x'' + 2D x' + (nu0^2 + nuc^2 cos(psi)) x = 0  for
%  (a) nu0^2 = nuc^2 = nu  (mean stiffness grows with the sweep), and
%  (b) nu0^2 = 1 fixed, nuc^2 = nu  (mean stiffness fixed, as for the
%      flapwise ODE, whose mean coefficients do not depend on mu).
clear; clc;

T = 2*pi; D = 0.15; N = 1024; m_range = -4:4;
nuv = linspace(0.05, 9, 46);
tf = linspace(0, T, N+1); tf(end) = [];
fr = -N/2:N/2-1;

cases = {'nu0^2 = nuc^2 = nu',        @(nu) nu, @(nu) nu;
         'nu0^2 = 1 fixed, nuc^2 = nu', @(nu) 1,  @(nu) nu};

for c = 1:size(cases,1)
    wmax = zeros(size(nuv));
    mdom = zeros(size(nuv));
    for k = 1:numel(nuv)
        a0 = cases{c,2}(nuv(k));
        ac = cases{c,3}(nuv(k));
        f = @(t,x) reshape([0 1; -(a0 + ac*cos(t)) -2*D]*reshape(x,2,2), 4, 1);
        sol = ode45(f, [0 T], reshape(eye(2),4,1), odeset('RelTol',1e-9,'AbsTol',1e-11));

        % Floquet exponent with positive imaginary part and its eigenvector
        PhiT = reshape(deval(sol,T),2,2);
        [V,L] = eig(PhiT);
        eta = log(diag(L))/T;
        [~,i] = max(imag(eta));

        % Periodic part of the mode (displacement) and its harmonics
        Ph = deval(sol, tf);
        Q = zeros(N,1);
        for j = 1:N
            Q(j) = [1 0]*reshape(Ph(:,j),2,2)*V(:,i)*exp(-eta(i)*tf(j));
        end
        C = fftshift(fft(Q)/N);
        p = arrayfun(@(m) abs(C(fr==m)), m_range);
        [~,im] = max(p);
        mdom(k) = m_range(im);
        wmax(k) = abs(m_range(im) + imag(eta(i)));
    end
    fprintf('%s\n  nu         : %s\n  k*         : %s\n  w(phi_max) : %s\n', cases{c,1}, ...
        mat2str(nuv(1:5:end),2), mat2str(mdom(1:5:end)), mat2str(wmax(1:5:end),2));
end
