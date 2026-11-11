% Cuffless Blood-Pressure via Pulse Transit Time (PTT) on real data
% =================================================================
% Estimates beat-to-beat systolic blood pressure from ECG + PPG alone and grades
% it against a real invasive arterial blood-pressure line.
%
% Data: BIDMC record 47 (PhysioNet), a real 8-minute ICU recording at 125 Hz with
% a simultaneous ECG lead (II), finger PPG (PLETH) and invasive arterial pressure
% (ABP). The ABP line is the clinical ground truth and is never used as an input.
%
% Chain: detect ECG R-peaks -> detect the PPG maximum-upslope point -> pulse
%        transit time per beat -> read true systolic pressure from ABP -> fit
%        SBP = a/PTT^2 + b on the first 90 s -> predict SBP for the rest with no
%        cuff -> report held-out error.
%
% The PTT reference point is the point of maximum upslope on the PPG, not the
% waveform foot. Both conventions are standard; maximum upslope is more robust to
% the baseline wander that dominates a finger PPG, and the two differ by a
% near-constant offset that the calibration absorbs.
%
% Both the R-peak and the upslope point are refined to sub-sample precision by a
% parabolic fit, because one sample at 125 Hz is 8 ms and a sample-indexed PTT
% would otherwise take only a handful of distinct values. Set USE_INTERP = false
% to reproduce the quantised version.
%
% This reproduces the Python implementation in cuffless_bp_ptt.py to the same
% numbers (450 beats, held-out MAE 3.0 mmHg, r = 0.45).
%
% Run in MATLAB (data/bidmc_47_Signals.csv must be alongside this file):
%   cuffless_bp_ptt
%
% Author: Yasar Aykut Altuntop

clear; close all;
fs = 125;
USE_INTERP = true;

T = readtable(fullfile('data','bidmc_47_Signals.csv'));
ecg = T.II;  ppg = T.PLETH;  abp = T.ABP;  t = (0:numel(ecg)-1)'/fs;

% R-peaks on the ECG
[bb,aa] = butter(3,[5 20]/(fs/2),'bandpass'); ecgf = filtfilt(bb,aa,ecg);
[~,rloc] = findpeaks(ecgf,'MinPeakDistance',round(0.4*fs),'MinPeakHeight',2.5*std(ecgf));

% PPG maximum-upslope point, on a low-passed PPG.
% gradient() is the central difference, matching numpy.gradient in the Python version.
[bl,al] = butter(4,8/(fs/2),'low'); ppgf = filtfilt(bl,al,ppg); d = gradient(ppgf);

ptt=[]; sbp=[]; tb=[];
for k = 1:numel(rloc)-1
    s0 = rloc(k); s1 = min(rloc(k)+round(0.5*fs), numel(t));
    if s1-s0 < 5, continue; end
    [~,rel] = max(d(s0:s1)); j = s0+rel-1;
    if USE_INTERP
        rAt = s0 + subsample_offset(ecgf, s0);
        pAt = j  + subsample_offset(d, j);
    else
        rAt = s0; pAt = j;
    end
    p = (pAt-rAt)/fs*1000;
    seg = abp(rloc(k):min(rloc(k)+round(fs),numel(abp)));
    sysP = max(seg);
    if p>50 && p<400 && sysP>60 && sysP<220
        ptt(end+1)=p; sbp(end+1)=sysP; tb(end+1)=t(rloc(k)); %#ok<SAGROW>
    end
end

% MAD outlier rejection + light median smoothing of the reference
madf = @(v) median(abs(v-median(v)))+1e-9;
keep = abs(ptt-median(ptt))<4*madf(ptt) & abs(sbp-median(sbp))<4*madf(sbp);
ptt=ptt(keep); sbp=sbp(keep); tb=tb(keep); sbp = medfilt1(sbp,5);

% calibrate on first 90 s, test on the rest
cal = tb <= tb(1)+90;
c   = polyfit(1./ptt(cal).^2, sbp(cal), 1);
est = polyval(c, 1./ptt.^2);
mae = mean(abs(est(~cal)-sbp(~cal)));
r   = corr(est(~cal)', sbp(~cal)');
fprintf('beats %d (cal %d / test %d)\n', numel(ptt), sum(cal), sum(~cal));
fprintf('mean PTT %.0f ms\n', mean(ptt));
fprintf('held-out MAE %.1f mmHg, r %.2f\n', mae, r);
fprintf('distinct PTT values: %d\n', numel(unique(round(ptt,6))));

figure; plot(tb, sbp, tb, est, '--', 'LineWidth',1.2);
xlabel('time [s]'); ylabel('SBP [mmHg]'); grid on;
legend('arterial line (truth)','cuffless PTT estimate');
title(sprintf('Cuffless SBP vs arterial line - held-out MAE %.1f mmHg, r %.2f', mae, r));


function off = subsample_offset(y, i)
%SUBSAMPLE_OFFSET  Sub-sample position of a local maximum at index i.
%   Parabolic fit through y(i-1), y(i), y(i+1). Returns the offset in samples
%   (|off| <= 1) so the refined location is i + off. Returns 0 at the array
%   edges or on a degenerate fit.
    off = 0;
    if i <= 1 || i >= numel(y), return; end
    denom = y(i-1) - 2*y(i) + y(i+1);
    if denom == 0, return; end
    o = 0.5*(y(i-1) - y(i+1))/denom;
    if abs(o) <= 1, off = o; end
end
