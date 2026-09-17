clc; clear; %close all;

load('participationTable.mat')

participation_data = participationTable{:,2:end};
nu_vals = participationTable{:,1};


mAll = size(participation_data,2);
mMax = (mAll-1)/2;


distPart = [(mMax-0.5):-1:0.5,0:mMax];
n = length(nu_vals);

iMax = zeros(n,1);
for k = 1:n
    [~, i0] = max(participation_data(k,:));
    cand = find(abs(participation_data(k,:) - participation_data(k,i0)) < 0.0025);
    [~, jj] = max(distPart(cand));
    iMax(k) = cand(jj);
end

participation_data_sel = participation_data(sub2ind(size(participation_data), (1:n)', iMax));

figure;
nexttile
plot(nu_vals, participation_data)
hold on; plot(nu_vals, participation_data_sel,'k:','linewidth',2)
%figure; plot(nu_vals,iMax)

domCurve = nan(length(nu_vals),mAll);
mask = zeros(length(nu_vals),mAll);
for idx = 1:mAll
    mask(:,idx) = (iMax == idx);
    domCurve(mask(:,idx) == 1, idx) = participation_data(find(mask(:,idx)== 1));
end

domCurve(mask, idx) = participation_data(mask);

plot(nu_vals, domCurve,'r*')

dominant_participation = max(participation_data');

nexttile
plot(nu_vals,distPart(iMax))
tmp = hampel(distPart(iMax));
hold on; plot(nu_vals,tmp,'r*');


