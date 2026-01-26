clc; clear; close all;
load('CharExTableRaw.mat','CharExTable');

Im_AllLine = CharExTable.Freq_s1_norm + (0:0.5:2.5).*ones(size(CharExTable.Freq_s1_norm));
omega      = sqrt(CharExTable.nu02);

figure; hold on;
plot(Im_AllLine);
plot(omega,'k','LineWidth',1.5);
hold off;

Im_m   = NaN(size(omega));
Im_m_d = NaN(size(omega));
idxmin = NaN(size(omega));

for idx = 1:length(omega)
    diff = Im_AllLine(idx,:) - omega(idx);     % no abs
    mask = diff <= 0;                          % only from below

    if idx == 1
        % first point: just take closest-from-below
        if any(mask)
            [Im_m(idx), idx_loc] = max(Im_AllLine(idx,mask));
            idx_all              = find(mask);
            idxmin(idx)          = idx_all(idx_loc);
            Im_m_d(idx)          = omega(idx) - Im_m(idx);
        else
            [Im_m_d(idx), idxmin(idx)] = min(Im_AllLine(idx,:) - omega(idx));
            Im_m(idx) = Im_AllLine(idx,idxmin(idx));
        end
    else
        % later points: from below AND not smaller than previous Im_m
        mask = mask & (Im_AllLine(idx,:) >= Im_m(idx-1));

        if any(mask)
            [Im_m(idx), idx_loc] = max(Im_AllLine(idx,mask));  % closest from below
            idx_all              = find(mask);
            idxmin(idx)          = idx_all(idx_loc);
            Im_m_d(idx)          = omega(idx) - Im_m(idx);
        else
            % if nothing satisfies monotonic + below, keep previous value
            Im_m(idx)   = Im_m(idx-1);
            Im_m_d(idx) = omega(idx) - Im_m(idx);
            idxmin(idx) = idxmin(idx-1);
        end
    end
end
<<<<<<< HEAD
hold on;
plot(Im_m,'b','LineWidth',1.5);      % blue curve, non-decreasing
=======

figure; hold on;
plot(Im_m,'b','LineWidth',1.5);      % blue curve, non-decreasing
plot(omega,'k--','LineWidth',1.5);   % dashed omega
>>>>>>> 16342c24dd3d857b752736907883e637769b53ac
xlabel('Index');
ylabel('Value');
grid on;
