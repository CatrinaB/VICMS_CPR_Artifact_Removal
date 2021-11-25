
clear all

% --- Load the data

load asystole_and_cpr__learn.mat  
cpr_data = data;
load pure_vf__learn.mat
vf_data = data; 
clear data

%%
% 
i = 1; j = 1; TT = 300;  
% snr = 0; alpha = sqrt(1 / 10^(snr/10));
%opt.exct = [1]; opt.disp = 'off'; opt.maxiter = 100;
%
% extract the signals
v = vf_data(i).ecg;  v = v(:);
c = cpr_data(j).ecg; c = c(:);
u = cpr_data(j).x_signal; u = u(:);
%
vf_t   = vf_data(i).time; vf_t = vf_t(:);
vf_ts  = (vf_t(end)  - vf_t(1)) / length(v);  % VF  sampling time
cpr_t  = cpr_data(j).time; cpr_t = cpr_t(:);    
cpr_ts = (cpr_t(end) - cpr_t(1)) / length(c); % CPR sampling time
%
% resample

v_resampled = resample(v, 100, round(1 / vf_ts));
c_resampled = resample(c, 100, round(1 / cpr_ts));
u_resampled = resample(u, 100, round(1 / cpr_ts));
clear v
clear c
clear u

% vf_t_resampled = linspace(vf_t(1), vf_t(end), length(v_resampled));
% cpr_t_resampled = linspace(cpr_t(1), cpr_t(end), length(c_resampled));

% subplot(121)
% plot(vf_t)
% subplot(122)
% plot(vf_t_resampled)
% 
% subplot(121)
% plot(vf_t, v)
% 
% hold on
% plot(vf_t_resampled, v_resampled)
% hold off 
% 
% subplot(122)
% plot(cpr_t, c)
% 
% hold on
% plot(cpr_t_resampled, c_resampled)
% hold off
%%
% truncate the longer signal
vf_n  = length(v_resampled);
cpr_n = length(c_resampled);


if cpr_n > vf_n
  c_trunc = c_resampled(1:vf_n);
  v_trunc = v_resampled(:);
  u_trunc = u_resampled(1:vf_n);
  T = linspace(0, vf_t(end) - vf_t(1), vf_n); % time vector;
  % linspace = a vector of N equally distanced elements
  % linspace(start, end, no_elem)
else
  v_trunc = v_resampled(1:cpr_n);
  c_trunc = c_resampled(:);
  u_trunc = u_resampled(:);
  T = linspace(0, cpr_t(end) - cpr_t(1), cpr_n);
end

clear v_resampled
clear c_resampled
clear u_resampled

% subplot(221)
% plot(cpr_t_resampled, c_resampled)
% 
% subplot(223)
% plot(T, c_trunc)
% 
% subplot(222)
% plot(vf_t_resampled, v_resampled)
% 
% subplot(224)
% plot(T, v_trunc)
%%
% detrend
c_detr = detrend(c_trunc,'constant'); % normalizes(?) data -> mean value is 0
v_detr = detrend(v_trunc,'constant');
u_detr = detrend(u_trunc, 'constant');

clear c_trunc
clear v_trunc
clear u_trunc;

% subplot(121)
% plot(c_detr)
% 
% hold on
% plot(c_trunc)
% hold off
% 
% subplot(122)
% plot(v_detr)
% 
% hold on
% plot(v_trunc)
% hold off

%%
% normalize
v_norm = v_detr / norm(v_detr);
c_norm = c_detr / norm(c_detr);
u_norm = u_detr / norm(u_detr);

clear v_detr
clear c_detr
clear u_detr

% subplot(121)
% plot(c_detr)
% 
% hold on
% plot(c_norm)
% hold off
% 
% subplot(122)
% plot(v_detr)
% 
% hold on
% plot(v_norm)
% hold off

%%

% cs = alpha * c_norm;
y = v_norm + c_norm;

% --- Low-pass FIR filtering

% figure
% pwelch(cs), hold on, pwelch(v_norm)
% ax = axis;  
% plot([0.075 0.075 ax(3:4)],'--')
% set(gca,'fontsize',20), xlabel('x'), ylabel('y'), title('t')

%%
figure
[pcs, fcs] = pwelch(c_norm);
[pv, fv] = pwelch(v_norm);

subplot(211)
plot(fcs,10*log10(pcs))
hold on
plot(fv, 10*log10(pv))
hold off

subplot(212)
plot(fcs, pcs)
hold on
plot(fv, pv)
hold off

xlabel('Frequency (Hz)')
ylabel('PSD (dB/Hz)')
legend('cpr', 'vf')
%%
figure
pwelch(c_norm)
hold on
pwelch(v_norm)
h = get(gca, 'Children');
set(h(1), 'Color', 'r');
hold off

legend('CPR', 'VF')

%%
% design low-pass filter
cpr_filter = fir1(70, 0.082); 
vf_filter = fir1(100, 0.074, 'high');

% Hd = dfilt.dffir(vf_filter);

%
% show magnitude response
% 
% fvtool(vf_filter), 
% fvtool(cpr_filter)
% set(gca,'fontsize',20), xlabel('x'), ylabel('y'), title('t')
%
%%
% figure

vh = filtfilt(vf_filter, 1, v_norm);
% subplot(221)
% plot(v_norm)
% subplot(223)
% plot(vh)
% 
% yt = fft(v_norm, 512);
% ys = fftshift(abs(yt));
% w = -256:255;
% subplot(222)
% plot(w, ys)
% yt = fft(vh, 512);
% ys = fftshift(abs(yt));
% subplot(224)
% plot(w, ys)
% 
% figure

ch = filtfilt(cpr_filter, 1, c_norm);
% subplot(221)
% plot(cs)
% subplot(223)
% plot(ch)
% 
% yt = fft(cs, 512);
% ys = fftshift(abs(yt));
% w = -256:255;
% subplot(222)
% plot(w, ys)
% yt = fft(ch, 512);
% ys = fftshift(abs(yt));
% subplot(224)
% plot(w, ys)


% fvtool(cs), fvtool(ch), fvtool(v_norm), fvtool(vh)

% save('c_and_v_normalized.mat', 'c_norm', 'v_norm', '-v6');

%%

ych = filtfilt(cpr_filter, 1, y); %filtered cpr
yvh = filtfilt(vf_filter, 1, y); %filtered vf

vyc = y - ych;
cyv = y - yvh;

% snr_v = 10 * log10(norm(v_norm) ^ 2 / norm(vh - v_norm) ^ 2)
snr_v_filt = 10 * log10(norm(v_norm) ^ 2 / norm(yvh - v_norm) ^ 2)
snr_v_yc = 10 * log10(norm(v_norm) ^ 2 / norm(vyc - v_norm) ^ 2)
% 10 * log10(norm(vh) ^ 2 / norm(yvh - vh) ^ 2)
% snr_c = 10 * log10(norm(cs) ^ 2 / norm(ch - cs) ^ 2 )
snr_c_filt = 10 * log10(norm(c_norm) ^ 2 / norm(ych - c_norm) ^ 2 )
snr_c_yv = 10 * log10(norm(c_norm) ^ 2 / norm(cyv - c_norm) ^ 2 )
% 10 * log10(norm(ch) ^ 2 / norm(ych - ch) ^ 2 )

% figure
% subplot(311)
% plot(y, 'LineWidth', 1), title('y')
% subplot(312)
% plot(ych, 'LineWidth', 1), title('CPR (filtered from y)')
% subplot(313)
% plot(yvh, 'LineWidth', 1), title('VF (filtered from y)')

% figure
% plot(yvh, 'LineWidth', 1.5), title('Signal comparison - VF')
% hold on
% plot(vh, 'r--', 'LineWidth', 1.5)
% plot(v_norm, 'g:', 'LineWidth', 2)
% hold off
% legend('VF filtered from y signal', 'Filtered VF input signal', 'Non-filtered VF input signal')
% axis([0 500 -0.1 0.1]), set(gca, 'fontsize', 15), grid minor

figure
plot(yvh, 'LineWidth', 1.5), title('Signal comparison - VF')
hold on
plot(v_norm, 'r--', 'LineWidth', 2)
plot(vyc, 'LineWidth', 1.5)
hold off
legend('VF filtered from y signal', 'Non-filtered VF input signal', 'v = y - c')
axis([0 500 -0.1 0.1]), set(gca, 'fontsize', 15), grid minor
% 
figure
plot(ych, 'LineWidth', 1.5), title('Signal comparison - CPR')
hold on
plot(c_norm, 'r--', 'LineWidth', 1.5)  
plot(cyv, 'LineWidth', 1.5)
hold off
legend('CPR filtered from y signal', 'Non-filtered CPR input signal', 'c = y - v')
axis([0 500 -0.1 0.1]), set(gca, 'fontsize', 15), grid minor


% save('c_and_v_normalized.mat', 'c_norm', 'v_norm', '-v6');