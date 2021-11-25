
clear all

% --- Load the data

load asystole_and_cpr__learn.mat  
cpr_data = data;
load pure_vf__learn.mat
vf_data = data; 
clear data

 
TT = 300;  m = 1; 
all_snr_v = zeros(49, 7); all_snr_c = zeros(49, 7);
all_snr_vyc = zeros(49, 7); all_snr_cyv = zeros(49, 7);
snr_y = [-10 -5 -0 5 10];
snr_v_metrics = zeros(5, 3); snr_c_metrics = zeros(5, 3); snr_vyc_metrics = zeros(5, 3);
snr_cyv_metrics = zeros(5, 3);
%%
for k=1:5
    alpha = snr_y(k)
    aaa = sqrt(1 / 10 ^ (alpha / 10))
    m=1;
    for i=1:7
        for j=1:7
            % extract the signals
            v = vf_data(i).ecg;  v = v(:);
            c = cpr_data(j).ecg; c = c(:);
            
            vf_t   = vf_data(i).time; vf_t = vf_t(:);
            vf_ts  = (vf_t(end)  - vf_t(1)) / length(v);  % VF  sampling time
            cpr_t  = cpr_data(j).time; cpr_t = cpr_t(:);    
            cpr_ts = (cpr_t(end) - cpr_t(1)) / length(c); % CPR sampling time
            
            % resample
            
            v_resampled = resample(v, 100, round(1 / vf_ts));
            c_resampled = resample(c, 100, round(1 / cpr_ts));
            clear v
            clear c
            
            % truncate the longer signal
            vf_n  = length(v_resampled);
            cpr_n = length(c_resampled);
            
            
            if cpr_n > vf_n
              c_trunc = c_resampled(1:vf_n);
              v_trunc = v_resampled(:);
              T = linspace(0, vf_t(end) - vf_t(1), vf_n); % time vector;
              % linspace = a vector of N equally distanced elements
              % linspace(start, end, no_elem)
            else
              v_trunc = v_resampled(1:cpr_n);
              c_trunc = c_resampled(:);
              T = linspace(0, cpr_t(end) - cpr_t(1), cpr_n);
            end
            
            clear v_resampled
            clear c_resampled
            
            % detrend
            c_detr = detrend(c_trunc,'constant'); % normalizes(?) data -> mean value is 0
            v_detr = detrend(v_trunc,'constant');
            clear c_trunc
            clear v_trunc
            
            % normalize
            v_norm = v_detr / norm(v_detr);
            c_norm = c_detr / norm(c_detr);
            
            clear v_detr
            clear c_detr
            
            y = v_norm + c_norm * aaa;
            
            % --- Low-pass FIR filtering
            
            % figure
            % pwelch(cs), hold on, pwelch(v_norm)
            % ax = axis;  
            % plot([0.075 0.075 ax(3:4)],'--')
            % set(gca,'fontsize',20), xlabel('x'), ylabel('y'), title('t')
                    
            % design low-pass filter
            cpr_filter = fir1(70, 0.082); 
            vf_filter = fir1(100, 0.074, 'high');
            
            %
            % show magnitude response
            % 
            % fvtool(vf_filter), 
            % fvtool(cpr_filter)
            
            ych = filtfilt(cpr_filter, 1, y); %filtered cpr
            yvh = filtfilt(vf_filter, 1, y); %filtered vf
            
            vyc = y - ych;
            cyv = y - yvh;
            
            all_snr_v(m, k) = 10 * log10(norm(v_norm) ^ 2 / norm(yvh - v_norm) ^ 2);
            all_snr_v(m, 6) = i; all_snr_v(m, 7) = j;
            all_snr_vyc(m, k) = 10 * log10(norm(v_norm) ^ 2 / norm(vyc - v_norm) ^ 2);
            all_snr_vyc(m, 6) = i; all_snr_vyc(m, 7) = j;
            
            all_snr_c(m, k) = 10 * log10(norm(c_norm) ^ 2 / norm(ych - c_norm) ^ 2 );
            all_snr_c(m, 6) = i; all_snr_c(m, 7) = j;
            all_snr_cyv(m, k) = 10 * log10(norm(c_norm) ^ 2 / norm(cyv - c_norm) ^ 2 );
            all_snr_cyv(m, 6) = i; all_snr_cyv(m, 7) = j;
    
            m = m + 1;
        end
    end
end

disp('end')
%%
for k=1:5
    snr_v_metrics(k, 1) = mean(all_snr_v(:, k));
    snr_v_metrics(k, 2) = max(all_snr_v(:, k));
    snr_v_metrics(k, 3) = min(all_snr_v(:, k));

    snr_c_metrics(k, 1) = mean(all_snr_c(:, k));
    snr_c_metrics(k, 2) = max(all_snr_c(:, k));
    snr_c_metrics(k, 3) = min(all_snr_c(:, k));

    snr_vyc_metrics(k, 1) = mean(all_snr_vyc(:, k));
    snr_vyc_metrics(k, 2) = max(all_snr_vyc(:, k));
    snr_vyc_metrics(k, 3) = min(all_snr_vyc(:, k));

    snr_cyv_metrics(k, 1) = mean(all_snr_cyv(:, k));
    snr_cyv_metrics(k, 2) = max(all_snr_cyv(:, k));
    snr_cyv_metrics(k, 3) = min(all_snr_cyv(:, k));

% maxv = max(all_snr_v(:,1));
% z = find(all_snr_v(:,1)==maxv, 1);
% fprintf('max snr v = %d, i = %d, j = %d \n', maxv, all_snr_v(z, 2), all_snr_v(z, 3))
% minv = min(all_snr_v(:,1));
% y = find(all_snr_v(:,1)==minv, 1)
% 
% fprintf('min snr v = %d, i = %d, j = %d \n', minv, all_snr_v(y, 2), all_snr_v(y, 3))


% avg_snr_c = mean(all_snr_c(:, 1));
% maxc = max(all_snr_c(:,1));
% z = find(all_snr_c(:,1)==maxc, 1)
% 
% fprintf('max snr c = %d, i = %d, j = %d \n', maxc, all_snr_c(z, 2), all_snr_c(z, 3))
% minc = min(all_snr_c(:,1));
% y = find(all_snr_c(:,1)==minc, 1)
% 
% fprintf('min snr c = %d, i = %d, j = %d \n', minc, all_snr_c(y, 2), all_snr_c(y, 3))
% 
% 
% avg_snr_vyc = mean(all_snr_vyc(:, 1));
% maxvyc = max(all_snr_vyc(:,1));
% z = find(all_snr_vyc(:,1)==maxvyc, 1)
% 
% fprintf('max snr vyc = %d, i = %d, j = %d \n', maxvyc, all_snr_vyc(z, 2), all_snr_vyc(z, 3))
% minvyc = min(all_snr_vyc(:,1));
% y = find(all_snr_vyc(:,1)==minvyc, 1)
% 
% fprintf('min snr vyc = %d, i = %d, j = %d \n', minvyc, all_snr_vyc(y, 2), all_snr_vyc(y, 3))
% 
% 
% avg_snr_cyv = mean(all_snr_cyv(:, 1));
% maxcyv = max(all_snr_cyv(:,1));
% z = find(all_snr_cyv(:,1)==maxcyv, 1)
% 
% fprintf('max snr cyv = %d, i = %d, j = %d \n', maxcyv, all_snr_cyv(z, 2), all_snr_cyv(z, 3))
% mincyv = min(all_snr_cyv(:,1));
% y = find(all_snr_cyv(:,1)==mincyv, 1)
% 
% fprintf('min snr cyv = %d, i = %d, j = %d \n', mincyv, all_snr_cyv(y, 2), all_snr_cyv(y, 3))
end
%%

% save('snr_c.mat', 'snr_c_metrics', '-v6');

% boxplot(all_snr_c(:, 1:5))
%% --- Kalman filtering     


