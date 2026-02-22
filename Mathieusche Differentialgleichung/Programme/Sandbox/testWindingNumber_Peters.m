clc; clear; close all;

%% Load and assign parameters
load('participation_data_and_M','participation_data', 'm_bubble','frequency_imag', 'nu_vals', 'growth_rate','composite_freq')
% omega = CharEx(:,1) .^ 0.5;
% nu2 = CharEx(:,1);
% vec = 0: 0.5:2.5;
% im_matrix = CharEx(:,5) + vec .* ones(size(CharEx(:,5)));

%% Calculate characteristic value closest to 

% Get the peaks of the composite_freq which is calculated from the frequncy
% branches and their modal participlation
[pksC,locsC] = findpeaks(composite_freq);
idx_m = [false; diff(m_bubble)~= 0];

% [sorted,idxSort] = sort(participation_data','descend');
% sorted2 = sorted(:,1:2);

% Look at the two largest branches to aoid jumps if two branches are nearly
% equal
sortedIndex = zeros(length(participation_data),2);
sortedIndexI = zeros(length(participation_data),2);
for idx = 1: length(participation_data)

    [sorted,idxSort] = sort(participation_data(idx,:),'descend');

    isNearlyEqual = sorted(1) - sorted(2) < 0.01;

   % if isNearlyEqual
        sortedIndex(idx,:) = sorted(1:2);
        sortedIndexI(idx,:) = sort(idxSort(1:2));
    % else
    %     sortedIndex(idx,1) = sorted(1);
    %     sortedIndexI(idx,1) = idxSort(1);
    % end
    % 

end

% Get the differenc between these sorted branches ( 1 by 500 vector)
diffSorted = diff(sortedIndex');

%Get the difference of that to be able to find the peaks
diffDiffSorted = [0,diff(diffSorted)];
[pksD,locsD] = findpeaks(diffDiffSorted,'MinPeakProminence',0.05);


% Combine information from the frequncy peak and the modal participation
tmp = zeros(size(m_bubble));
% First shift winding number from 
tmp(locsD(1)) = 0.5; 
tmp(locsC) = 0.5;
m_modpart = cumsum(tmp);


% Calculate the smoothed Peters winding number
[maxval, idxval ]= max(participation_data'); %#ok<UDIM>
[pks,locs] = findpeaks(maxval);
[pksNeg,locsNeg] = findpeaks(-maxval);
m_range = -4:4;
strLeg = arrayfun(@(m) sprintf('m=%d', m_range(m)), 1:length(m_range), 'UniformOutput', false);

%egend(arrayfun(@(m) sprintf('m=%d', m_range(m)), 1:length(m_range), 'UniformOutput', false), 'Location','northeastoutside');
cl = lines;
figure;
subplot(4,1,1)
plot(nu_vals,composite_freq);
hold on; plot(nu_vals(locsC),pksC,'k*');
plot(nu_vals(idx_m),composite_freq(idx_m),'ro');
plot(nu_vals(locsD),composite_freq(locsD),'bo','MarkerSize',14);
axis tight;
grid on;

ylabel('Frequency \omega')
subplot(4,1,2)
plot(nu_vals, participation_data);
hold on; plot(nu_vals, maxval,'r--','LineWidth',1.3)

%plot(locs,pks,'*')
%plot(nu_vals(locsNeg),-pksNeg,'ko')
ylabel('mod part.')
axis tight;
grid on;

subplot(4,1,3)
plot(nu_vals,sortedIndex(:,1),'color',cl(1,:));
hold on;
plot(nu_vals,sortedIndex(:,2),'--','LineWidth',1.3,'color',cl(2,:));
plot(nu_vals,diffDiffSorted,'k');
grid on;
axis tight;

subplot(4,1,4)
plot(nu_vals,m_bubble,'color',cl(1,:));
hold on;
plot(nu_vals, m_modpart,'--','color',cl(2,:))
ylabel('m bubble counting')


axis tight;
grid on;



xlabel('Amplification factor \nu')


%% Unused code, keep for later
% figure;
% plot(nu_vals,idxval)
% ylabel('max value participation')

% Look at the two largest branches to aoid jumps if two branches are nearly
% % equal
% sortedIndex = zeros(length(participation_data),2);
% sortedIndexI = zeros(length(participation_data),2);
% for idx = 1: length(participation_data)
% 
%     [sorted,idxSort] = sort(participation_data(idx,:),'descend');
% 
%     isNearlyEqual = sorted(1) - sorted(2) < 0.01;
% 
%    % if isNearlyEqual
%         sortedIndex(idx,:) = sorted(1:2);
%         sortedIndexI(idx,:) = sort(idxSort(1:2));
%     % else
%     %     sortedIndex(idx,1) = sorted(1);
%     %     sortedIndexI(idx,1) = idxSort(1);
%     % end
%     % 

end

