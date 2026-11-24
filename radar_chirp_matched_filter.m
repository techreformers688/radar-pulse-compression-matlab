%% ========================================================================
%  RADAR CHIRP SIMULATION WITH MATCHED FILTERING (PULSE COMPRESSION)
%  Author : Srishti Tiwari
%  Tool   : MATLAB. Core MATLAB only; the Phased Array System Toolbox is
%           used ONLY for optional cross-checks (wrapped in try/catch).
%  Run    : press Run. Figures are saved to ./figures, numbers to results.txt
%
%  What this script demonstrates (each section is one claim you can defend):
%   1  Radar parameters and theory values to compare against
%   2  LFM chirp generated from first principles (no toolbox)
%   3  Noisy multi-target echo where the raw pulse is invisible
%   4  Matched filter in time domain AND frequency domain (they must agree)
%   5  Measured SNR gain, range resolution and sidelobes vs theory
%   6  Resolution vs bandwidth (two close targets merge, then split)
%   7  Sidelobe suppression with a window, and what it costs
%   8  Doppler: range-Doppler coupling and ambiguity function
%   9  Detection probability vs SNR: raw threshold vs matched filter
%  10  Cross-check against the Phased Array System Toolbox (optional)
%  11  Results table (theory vs measured), saved to results.txt
%% ========================================================================
clear; clc; close all; rng(42);
scriptDir = fileparts(mfilename('fullpath'));
figDir = fullfile(scriptDir, 'figures');
if ~exist(figDir, 'dir'), mkdir(figDir); end


%% 1. RADAR PARAMETERS ----------------------------------------------------
c   = 3e8;                 % speed of light (m/s)
B   = 10e6;                % LFM sweep bandwidth (Hz)
T   = 20e-6;               % pulse width (s)
PRI = 100e-6;              % pulse repetition interval (s)
OS  = 8;                   % oversampling factor
fs  = OS*B;                % complex baseband sample rate (Hz)
k   = B/T;                 % chirp rate (Hz/s)
Ns  = round(T*fs);         % samples per pulse
N   = round(PRI*fs);       % samples per PRI
TB  = B*T;                 % time-bandwidth product

snr_in_dB = -10;           % input SNR measured in the signal bandwidth B
% Signal power is 1, so noise variance per sample (noise spread over fs):
sigma2 = (fs/B) / 10^(snr_in_dB/10);

% Theory values to compare against
th.R_unamb = c*PRI/2;                 % unambiguous range (m)
th.res_nom = c/(2*B);                 % nominal range resolution (m)
th.res_3dB = 0.886*c/(2*B);           % -3 dB mainlobe width, unweighted LFM (m)
th.gain_dB = 10*log10(TB);            % pulse compression gain (dB)
th.psl_dB  = -13.26;                  % sinc peak sidelobe level (dB)

fprintf('Chirp rate k = %.3g Hz/s, TB = %d, unambiguous range = %.1f km\n', ...
        k, TB, th.R_unamb/1e3);

%% 2. GENERATE THE LFM CHIRP (from first principles) -----------------------
[s, tt] = make_chirp(B, T, fs);        % complex baseband, sweeps -B/2 ... +B/2

f_inst = diff(unwrap(angle(s))) * fs/(2*pi);       % instantaneous frequency
NF = 8192;
S  = fftshift(fft(s, NF));
f  = (-NF/2:NF/2-1)' * fs/NF;
Sdb = 20*log10(abs(S)/max(abs(S)));
bw_meas = max(f(Sdb > -10)) - min(f(Sdb > -10));   % -10 dB occupied bandwidth

figure('Name','LFM chirp','Position',[100 100 1000 700]);
subplot(3,1,1); plot(tt*1e6, real(s)); grid on;
xlabel('Time (\mus)'); ylabel('Re\{s(t)\}'); title('LFM chirp, real part');
subplot(3,1,2); plot(tt(1:end-1)*1e6, f_inst/1e6); grid on;
xlabel('Time (\mus)'); ylabel('Frequency (MHz)');
title(sprintf('Instantaneous frequency sweeps %.0f MHz (k = %.2g Hz/s)', B/1e6, k));
subplot(3,1,3); plot(f/1e6, Sdb); grid on; ylim([-40 5]); xlim([-fs/2e6/2 fs/2e6/2]);
xlabel('Frequency (MHz)'); ylabel('|S(f)| (dB)');
title(sprintf('Spectrum: measured -10 dB bandwidth = %.2f MHz (design B = %.0f MHz)', ...
      bw_meas/1e6, B/1e6));
saveas(gcf, fullfile(figDir,'01_chirp.png'));

%% 3. NOISY MULTI-TARGET ECHO ---------------------------------------------
ranges = [3000 5000 5010];            % m  (last two are 10 m apart)
amps   = [1.0  0.8  0.8];
rx_clean = zeros(N,1);
for i = 1:numel(ranges)
    rx_clean = rx_clean + amps(i)*make_echo(2*ranges(i)/c, fs, N, B, T, 0);
end
rx = rx_clean + noise_matrix(sigma2, N, 1);
r  = (0:N-1)' / fs * c/2;             % range axis (m)

%% 4. MATCHED FILTER: TIME DOMAIN vs FREQUENCY DOMAIN ----------------------
rect = ones(Ns,1);
y_f = mfilt(rx, s, rect);                           % frequency-domain (FFT) implementation
y_t = conv(rx, conj(flipud(s)));                    % time-domain implementation
y_t = y_t(Ns:Ns+N-1);                               % align lag 0
impl_err = max(abs(y_f - y_t)) / max(abs(y_f));
fprintf('Time vs frequency domain matched filter, max relative error = %.2e\n', impl_err);

figure('Name','Pulse compression','Position',[100 100 1000 700]);
subplot(2,1,1); plot(r/1e3, 20*log10(abs(rx)+eps)); grid on;
xlabel('Range (km)'); ylabel('dB'); xlim([0 th.R_unamb/1e3]);
title(sprintf('Raw received signal, input SNR = %d dB: targets are not visible', snr_in_dB));
subplot(2,1,2); plot(r/1e3, 20*log10(abs(y_f)+eps)); grid on; hold on;
yl = ylim; for i = 1:numel(ranges), plot(ranges(i)/1e3*[1 1], yl, 'r--'); end
xlabel('Range (km)'); ylabel('dB'); xlim([0 th.R_unamb/1e3]);
title('After matched filter: targets at 3.0 km and ~5.0 km appear (red = truth)');
saveas(gcf, fullfile(figDir,'02_raw_vs_compressed.png'));

%% 5. MEASURED SNR GAIN, RESOLUTION, SIDELOBES vs THEORY -------------------
x1    = make_echo(2*5000/c, fs, N, B, T, 0);        % single noiseless target at 5 km
wins  = {rect, hamm(Ns)};
wname = {'Rectangular (no window)', 'Hamming'};
meas  = struct();
for i = 1:2
    w  = wins{i};
    y1 = mfilt(x1, s, w);
    [pk, ipk] = max(abs(y1));
    yn = mfilt(noise_matrix(sigma2, N, 100), s, w);       % noise-only output, 100 runs
    nz = mean(mean(abs(yn(1:N-Ns,:)).^2));                % exclude truncated tail
    meas(i).snr_out_dB = 10*log10(pk^2/nz);
    meas(i).gain_dB    = meas(i).snr_out_dB - snr_in_dB;
    [meas(i).width_m, meas(i).psl_dB] = pulse_metrics(y1, fs, c);
    meas(i).peak_err_m = (ipk-1)/fs*c/2 - 5000;
    meas(i).y = y1;
end
fprintf('\n%-26s gain %.2f dB | -3dB width %.2f m | PSL %.1f dB\n', wname{1}, ...
        meas(1).gain_dB, meas(1).width_m, meas(1).psl_dB);
fprintf('%-26s gain %.2f dB | -3dB width %.2f m | PSL %.1f dB\n', wname{2}, ...
        meas(2).gain_dB, meas(2).width_m, meas(2).psl_dB);
snr_loss_dB = meas(1).gain_dB - meas(2).gain_dB;

%% 6. RESOLUTION vs BANDWIDTH ---------------------------------------------
[rA, aA] = pair_response(10e6, T, OS, PRI, c, 5000, 5010);
[rB, aB] = pair_response(30e6, T, OS, PRI, c, 5000, 5010);
nA = count_peaks(rA, aA, 4950, 5060);
nB = count_peaks(rB, aB, 4950, 5060);

Bs = (2:2:40)' * 1e6;
wid = zeros(size(Bs));
for i = 1:numel(Bs)
    fsx = OS*Bs(i);
    sx  = make_chirp(Bs(i), T, fsx);
    z   = zeros(numel(sx),1);
    yx  = mfilt([z; sx; z], sx, ones(numel(sx),1));
    wid(i) = pulse_metrics(yx, fsx, c);
end

figure('Name','Resolution','Position',[100 100 1000 700]);
subplot(2,1,1); hold on;
plot(rA, 20*log10(aA/max(aA)+eps), 'LineWidth',1.5);
plot(rB, 20*log10(aB/max(aB)+eps), 'LineWidth',1.5);
xlim([4950 5060]); ylim([-30 3]); grid on;
for rr = [5000 5010], plot(rr*[1 1], [-30 3], 'k:'); end
legend(sprintf('B = 10 MHz (%.0f m res): %d peak', c/(2*10e6), nA), ...
       sprintf('B = 30 MHz (%.0f m res): %d peaks', c/(2*30e6), nB), 'Location','south');
xlabel('Range (m)'); ylabel('Normalised (dB)');
title('Two equal targets 10 m apart: unresolved at 10 MHz, resolved at 30 MHz');
subplot(2,1,2); hold on;
plot(Bs/1e6, wid, 'o-', 'LineWidth',1.5);
plot(Bs/1e6, 0.886*c./(2*Bs), '--', 'LineWidth',1.5);
plot(Bs/1e6, c./(2*Bs), ':', 'LineWidth',1.5); grid on;
legend('Measured -3 dB width','Theory 0.886 c/2B','Nominal c/2B','Location','northeast');
xlabel('Bandwidth (MHz)'); ylabel('Range resolution (m)');
title('Range resolution vs bandwidth (independent of pulse width)');
saveas(gcf, fullfile(figDir,'03_resolution.png'));

%% 7. SIDELOBES: WINDOWING AND ITS COST ------------------------------------
x_strong = make_echo(2*4000/c, fs, N, B, T, 0);
x_weak   = 10^(-25/20) * make_echo(2*4060/c, fs, N, B, T, 0);   % -25 dB, 60 m away
x_two    = x_strong + x_weak;
y_two_r  = mfilt(x_two, s, wins{1});
y_two_h  = mfilt(x_two, s, wins{2});

figure('Name','Sidelobes','Position',[100 100 1000 700]);
subplot(2,1,1); hold on;
y1r = meas(1).y; y1h = meas(2).y;
plot(r, 20*log10(abs(y1r)/max(abs(y1r))+eps), 'LineWidth',1.2);
plot(r, 20*log10(abs(y1h)/max(abs(y1h))+eps), 'LineWidth',1.2);
xlim([4850 5150]); ylim([-70 3]); grid on;
legend(sprintf('Rectangular: PSL %.1f dB, width %.1f m', meas(1).psl_dB, meas(1).width_m), ...
       sprintf('Hamming: PSL %.1f dB, width %.1f m', meas(2).psl_dB, meas(2).width_m));
xlabel('Range (m)'); ylabel('Normalised (dB)');
title(sprintf('Windowing trades resolution and %.2f dB of SNR for lower sidelobes', snr_loss_dB));
subplot(2,1,2); hold on;
plot(r, 20*log10(abs(y_two_r)/max(abs(y_two_r))+eps), 'LineWidth',1.2);
plot(r, 20*log10(abs(y_two_h)/max(abs(y_two_h))+eps), 'LineWidth',1.2);
xlim([3900 4200]); ylim([-70 3]); grid on; yl = ylim;
plot(4060*[1 1], yl, 'k:');
legend('Rectangular','Hamming','Weak target (-25 dB) at 4060 m');
xlabel('Range (m)'); ylabel('Normalised (dB)');
title('Strong target at 4000 m and a weak one at 4060 m: sidelobes vs masking');
saveas(gcf, fullfile(figDir,'04_sidelobes.png'));

%% 8. DOPPLER: RANGE-DOPPLER COUPLING AND AMBIGUITY FUNCTION ---------------
fd_list = [0 100e3 200e3];
dop = struct();
figure('Name','Doppler','Position',[100 100 1000 700]);
subplot(2,2,1); hold on;
for i = 1:numel(fd_list)
    xd = make_echo(2*5000/c, fs, N, B, T, fd_list(i));
    yd = mfilt(xd, s, rect);
    [pk, ip] = max(abs(yd));
    dop(i).shift_m  = (ip-1)/fs*c/2 - 5000;
    dop(i).theory_m = fd_list(i)/k * c/2;
    dop(i).loss_dB  = 20*log10(pk/Ns);
    plot(r, 20*log10(abs(yd)/Ns+eps), 'LineWidth',1.2);
end
xlim([4900 5100]); ylim([-40 3]); grid on;
legend(arrayfun(@(v) sprintf('f_d = %d kHz', v/1e3), fd_list, 'UniformOutput', false));
xlabel('Range (m)'); ylabel('dB'); title('Peak shifts with Doppler (range-Doppler coupling)');

fd_axis = (-100e3:2e3:100e3)';
lag     = (-(Ns-1):(Ns-1))' / fs;
hrev    = conj(flipud(s));
AF      = zeros(numel(lag), numel(fd_axis));
for i = 1:numel(fd_axis)
    AF(:,i) = abs(conv(s .* exp(1j*2*pi*fd_axis(i)*tt), hrev)) / Ns;
end
AFdb = max(20*log10(AF + 1e-6), -40);
subplot(2,2,2);
imagesc(fd_axis/1e3, lag*1e6, AFdb); axis xy; colorbar;
try, clim([-40 0]); catch, caxis([-40 0]); end
xlabel('Doppler (kHz)'); ylabel('Delay (\mus)'); title('Ambiguity function |\chi(\tau,f_d)| (dB)');
subplot(2,2,3);
[~, i0] = min(abs(fd_axis));            % column for zero Doppler
plot(lag*1e6, AFdb(:, i0), 'LineWidth',1.2); grid on; ylim([-40 3]);
xlabel('Delay (\mus)'); ylabel('dB'); title('Zero-Doppler cut (range profile)');
subplot(2,2,4);
plot(fd_axis/1e3, AFdb(Ns, :), 'LineWidth',1.2); grid on; ylim([-40 3]);
xlabel('Doppler (kHz)'); ylabel('dB'); title('Zero-delay cut (Doppler profile)');
saveas(gcf, fullfile(figDir,'05_doppler_ambiguity.png'));

% Optional cross-check of my zero-Doppler cut against ambgfun
amb_err = NaN;
try
    [afmag, delay_tb] = ambgfun(s, fs, 1/T, 'Cut', 'Delay');
    mine = interp1(lag, AF(:, i0), delay_tb, 'linear', 0);
    amb_err = max(abs(mine(:) - afmag(:)));
catch
end

%% 9. DETECTION: PD vs SNR, RAW THRESHOLD vs MATCHED FILTER ----------------
% SNR here is PER SAMPLE (noise variance 1, target amplitude A). The matched
% filter integrates all Ns samples coherently, so against a single-sample
% threshold the theoretical gap is 10*log10(Ns) = TB gain + oversampling (OS)
% effect. (TB gain, 23 dB, is quoted against SNR in the signal bandwidth B.)
Pfa = 1e-3;  M = 3000;
sg2 = 1;
thr_raw = -sg2 * log(Pfa);                % |x|^2 is exponential(sg2)
thr_mf  = -sg2 * Ns * log(Pfa);           % |s'x|^2 is exponential(sg2*Ns)
snr_sweep = -44:2:14;
gap_theory_dB = 10*log10(Ns);
Pd_raw = zeros(size(snr_sweep));  Pd_mf = Pd_raw;
Nz = sqrt(sg2/2) * (randn(Ns,M) + 1j*randn(Ns,M));
pfa_raw = mean(abs(Nz(round(Ns/2),:)).^2 > thr_raw);
pfa_mf  = mean(abs(s' * Nz).^2 > thr_mf);
for i = 1:numel(snr_sweep)
    A  = 10^(snr_sweep(i)/20);
    X  = Nz + A*s;                          % target present, same noise draws
    Pd_raw(i) = mean(abs(X(round(Ns/2),:)).^2 > thr_raw);
    Pd_mf(i)  = mean(abs(s' * X).^2 > thr_mf);
end
figure('Name','Detection','Position',[100 100 900 500]);
plot(snr_sweep, Pd_raw, 's-', snr_sweep, Pd_mf, 'o-', 'LineWidth',1.5); grid on;
legend('Naive threshold on raw samples','Matched filter + threshold','Location','southeast');
xlabel('Per-sample input SNR (dB)'); ylabel('P_d');
title(sprintf('Detection probability at P_{fa} = %.0e (%d Monte Carlo trials per point)', Pfa, M));
saveas(gcf, fullfile(figDir,'06_detection.png'));
snr_pd50_raw = interp_cross(snr_sweep, Pd_raw, 0.5);
snr_pd50_mf  = interp_cross(snr_sweep, Pd_mf, 0.5);

%% 10. CROSS-CHECK WITH THE PHASED ARRAY SYSTEM TOOLBOX (optional) ----------
w_tb = NaN; psl_tb = NaN; tb_ok = false;
try
    wfm = phased.LinearFMWaveform('SampleRate',fs, 'PulseWidth',T, 'PRF',1/PRI, ...
          'SweepBandwidth',B, 'SweepDirection','Up');
    xtb = wfm();
    mf  = phased.MatchedFilter('Coefficients', getMatchedFilter(wfm));
    ytb = mf(xtb);
    [w_tb, psl_tb] = pulse_metrics(ytb, fs, c);
    tb_ok = true;
catch
end

%% 11. RESULTS TABLE (theory vs measured) -----------------------------------
L = {};
L{end+1} = sprintf('RADAR CHIRP + MATCHED FILTER: RESULTS');
L{end+1} = sprintf('B = %.0f MHz, T = %.0f us, TB = %d, fs = %.0f MHz, PRI = %.0f us, input SNR = %d dB', ...
                   B/1e6, T*1e6, TB, fs/1e6, PRI*1e6, snr_in_dB);
L{end+1} = repmat('-', 1, 78);
L{end+1} = sprintf('%-44s %12s %12s', 'Metric', 'Theory', 'Measured');
L{end+1} = sprintf('%-44s %12.2f %12.2f', 'Pulse compression gain (dB)', th.gain_dB, meas(1).gain_dB);
L{end+1} = sprintf('%-44s %12.2f %12.2f', '-3 dB range resolution, rectangular (m)', th.res_3dB, meas(1).width_m);
L{end+1} = sprintf('%-44s %12.2f %12s',   'Nominal range resolution c/2B (m)', th.res_nom, '(reference)');
L{end+1} = sprintf('%-44s %12.1f %12.1f', 'Peak sidelobe, rectangular (dB)', th.psl_dB, meas(1).psl_dB);
L{end+1} = sprintf('%-44s %12s %12.1f', 'Peak sidelobe, Hamming (dB)', '-', meas(2).psl_dB);
L{end+1} = sprintf('%-44s %12s %12.2f', 'Hamming mainlobe width (m)', '-', meas(2).width_m);
L{end+1} = sprintf('%-44s %12s %12.2f', 'Hamming SNR loss vs rectangular (dB)', '-', snr_loss_dB);
L{end+1} = sprintf('%-44s %12s %12.2e', 'Time vs frequency-domain MF rel. error', '-', impl_err);
L{end+1} = sprintf('%-44s %12s %12.2f', 'Measured -10 dB bandwidth (MHz)', sprintf('%.2f', B/1e6), bw_meas/1e6);
L{end+1} = sprintf('%-44s %12s %12d', 'Peaks found, 10 m pair @ 10 MHz', '1', nA);
L{end+1} = sprintf('%-44s %12s %12d', 'Peaks found, 10 m pair @ 30 MHz', '2', nB);
for i = 1:numel(fd_list)
    L{end+1} = sprintf('%-44s %12.1f %12.1f', sprintf('Range shift at fd = %d kHz (m)', fd_list(i)/1e3), ...
                       abs(dop(i).theory_m), abs(dop(i).shift_m)); %#ok<AGROW>
end
L{end+1} = sprintf('%-44s %12s %12.2f', 'Peak loss at f_d = 200 kHz (dB)', '-', dop(3).loss_dB);
L{end+1} = sprintf('%-44s %12.0e %12.1e', 'Empirical Pfa, matched filter (target Pfa)', Pfa, pfa_mf);
L{end+1} = sprintf('%-44s %12s %12.1f', 'Per-sample SNR for Pd = 0.5, raw sample (dB)', '-', snr_pd50_raw);
L{end+1} = sprintf('%-44s %12s %12.1f', 'Per-sample SNR for Pd = 0.5, matched filter (dB)', '-', snr_pd50_mf);
L{end+1} = sprintf('%-44s %12.2f %12.1f', 'Detection gap raw vs MF at Pd = 0.5 (dB)', gap_theory_dB, snr_pd50_raw - snr_pd50_mf);
if tb_ok
    L{end+1} = sprintf('%-44s %12s %12.2f', 'Toolbox MatchedFilter: -3 dB width (m)', sprintf('%.2f', meas(1).width_m), w_tb);
    L{end+1} = sprintf('%-44s %12s %12.1f', 'Toolbox MatchedFilter: PSL (dB)', sprintf('%.1f', meas(1).psl_dB), psl_tb);
else
    L{end+1} = 'Toolbox cross-check skipped (Phased Array System Toolbox not available).';
end
if ~isnan(amb_err)
    L{end+1} = sprintf('%-44s %12s %12.2e', 'My ambiguity cut vs ambgfun, max abs diff', '-', amb_err);
end
fprintf('\n%s\n', strjoin(L, newline));
fid = fopen(fullfile(scriptDir,'results.txt'),'w'); fprintf(fid, '%s\n', strjoin(L, newline)); fclose(fid);
fprintf('\nSaved figures to %s and the table to %s\n', figDir, fullfile(scriptDir,'results.txt'));

%% ========================================================================
%  LOCAL FUNCTIONS
%% ========================================================================
function [s, tt] = make_chirp(B, T, fs)
% Complex baseband LFM chirp sweeping -B/2 ... +B/2 over T seconds.
    Ns = round(T*fs);
    tt = (0:Ns-1)' / fs;
    k  = B/T;
    s  = exp(1j*pi*k*(tt - T/2).^2);
end

function x = make_echo(tau, fs, N, B, T, fd)
% Received baseband echo of the chirp delayed by tau (fractional delays OK,
% no sample quantisation) with an optional Doppler shift fd.
    t  = (0:N-1)' / fs;
    u  = t - tau;
    on = (u >= 0) & (u < T);
    x  = zeros(N,1);
    x(on) = exp(1j*pi*(B/T)*(u(on) - T/2).^2) .* exp(1j*2*pi*fd*t(on));
end

function nz = noise_matrix(sigma2, N, M)
% Circular complex Gaussian noise, variance sigma2 per sample.
    nz = sqrt(sigma2/2) * (randn(N,M) + 1j*randn(N,M));
end

function y = mfilt(rx, s, w)
% Frequency-domain matched filter. Reference replica = chirp .* window w.
% Output index n (0-based) is the correlation at a delay of n samples.
    N = size(rx,1);  Ns = numel(s);
    Nfft = 2^nextpow2(N + Ns - 1);
    h = s(:) .* w(:);
    y = ifft(fft(rx, Nfft) .* conj(fft(h, Nfft)));
    y = y(1:N, :);
end

function w = hamm(n)
    w = 0.54 - 0.46*cos(2*pi*(0:n-1)'/(n-1));
end

function [width_m, psl_db] = pulse_metrics(y, fs, c)
% -3 dB mainlobe width (m) and peak sidelobe level (dB re peak) of a
% compressed pulse, with linear interpolation at the -3 dB crossings.
    a = abs(y(:));  [pk, ip] = max(a);
    ydb = 20*log10(a/pk + eps);  lvl = 20*log10(1/sqrt(2));
    il = ip;  while il > 1 && ydb(il) > lvl, il = il - 1; end
    xl = il + (lvl - ydb(il)) / (ydb(il+1) - ydb(il));
    ir = ip;  while ir < numel(ydb) && ydb(ir) > lvl, ir = ir + 1; end
    xr = ir - 1 + (ydb(ir-1) - lvl) / (ydb(ir-1) - ydb(ir));
    width_m = (xr - xl) / fs * c/2;
    nl = ip;  while nl > 1 && a(nl-1) < a(nl), nl = nl - 1; end      % first null, left
    nr = ip;  while nr < numel(a) && a(nr+1) < a(nr), nr = nr + 1; end % first null, right
    psl_db = max([ydb(1:nl-1); ydb(nr+1:end)]);
end

function [r, a] = pair_response(Bx, T, OS, PRI, c, R1, R2)
% Compressed response of two equal targets at ranges R1, R2 for bandwidth Bx.
    fsx = OS*Bx;  sx = make_chirp(Bx, T, fsx);  Nx = round(PRI*fsx);
    x = make_echo(2*R1/c, fsx, Nx, Bx, T, 0) + make_echo(2*R2/c, fsx, Nx, Bx, T, 0);
    y = mfilt(x, sx, ones(numel(sx),1));
    r = (0:Nx-1)' / fsx * c/2;   a = abs(y);
end

function n = count_peaks(r, a, rlo, rhi)
% Number of local maxima within [rlo,rhi] that are within 6 dB of the largest.
    idx = find(r >= rlo & r <= rhi);
    a = a / max(a(idx));  n = 0;
    for i = idx(2:end-1)'
        if a(i) > a(i-1) && a(i) >= a(i+1) && a(i) > 10^(-6/20), n = n + 1; end
    end
end

function xc = interp_cross(x, y, level)
% x where the (increasing) curve y first crosses 'level'.
    i = find(y >= level, 1, 'first');
    if isempty(i), xc = NaN; elseif i == 1, xc = x(1);
    else, xc = x(i-1) + (level - y(i-1)) * (x(i) - x(i-1)) / (y(i) - y(i-1)); end
end
