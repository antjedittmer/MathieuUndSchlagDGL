clc; clear; close all;

%% Load and assign parameters
load('CharExMatrix.mat','CharEx');
omega = CharEx(:,1) .^ 0.5;
nu2 = CharEx(:,1);
vec = 0: 0.5:2.5;
im_matrix = CharEx(:,5) + vec .* ones(size(CharEx(:,5)));

%% Calculate characteristic value closest to 
im_min = nan(length(omega),1);
im_idx  = nan(length(omega),1);
for idx = 1: length(omega)
    tmpIdx  = omega(idx) >= im_matrix(idx,:); % only use values smaller than omega
    tmp_im_matrix = im_matrix(idx,tmpIdx); % values from currem row smaller than omega
    [~,im_idx(idx)] = min(abs(omega(idx) - tmp_im_matrix)); % index value close to omega
    im_min(idx) = im_matrix(idx,im_idx(idx)); %imag. part char. exp. with addition factor closest to omega
end

figure;
tiledlayout(2,1,'TileSpacing','compact');

nexttile
plot(nu2,im_matrix); hold on; 
plot(nu2,omega,'k')

hold on;
plot(nu2,im_min,'r--','LineWidth',1.2);
axis tight; grid on;

ylabel('Imag. part characteristic expo.');
tmp = regexp(sprintf('Imag. part + %2.1f,', vec ),',','split');
legCell = [tmp(1:end-1),'\omega','Imag. part + m'];
legend(legCell,'Location','eastoutside');

nexttile
plot(nu2,vec(im_idx));
ylabel('Addition term m')
legend('Term m','Location','eastoutside')
xlabel('Amplification factor \nu_C^2 = \omega^2');
axis tight; grid on;

exportgraphics(gcf,'test_compare_imagCharExp_omega.png','Resolution',300)


