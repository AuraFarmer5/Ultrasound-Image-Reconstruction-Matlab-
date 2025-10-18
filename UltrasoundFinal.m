close all
clear all

%Import CSV File
rawData = readtable("TiO2 D4.csv");
clc

%% Extract relevant columns from the imported data
T = rawData;
time = double(T.WaveformType(:));
signal = double(T.ANALOG(:));
synch = double(T.Var3(:));

%% Clean & sort spiking
ok = isfinite(time) & isfinite(signal) & isfinite(synch);
time  = time(ok);  signal = signal(ok);  synch = synch(ok);
[time, idx] = sort(time);
signal = signal(idx);  synch = synch(idx);
fs = 1/mean(diff(time));

% Set t=0 at SYNC rising edge 
thr = (max(synch)+min(synch))/2;
k0  = find(synch(1:end-1)<thr & synch(2:end)>=thr, 1, 'first');
if ~isempty(k0)
    t0 = time(k0);
else
    warning('SYNC rising edge not found; using min(time) as zero.');
    t0 = min(time);
end
time = time - t0;

% Transient handling: short blank + optional fade 
blank_us = 3; % try 2–5 µs
Nblank   = round(blank_us*1e-6*fs);
i0 = max(k0,1); 
i1 = min(i0+Nblank-1, numel(signal));

signal2 = signal;
if ~isempty(k0)
    signal2(i0:i1) = 0; % mute transient
    % 0.5 µs fade-in
    Nfade = round(0.5e-6*fs);
    j0 = min(i1+1, numel(signal));
    j1 = min(i1+Nfade, numel(signal));
    if j1>j0
        w = hann(2*Nfade); w = w(Nfade+1:end);
        signal2(j0:j1) = signal(j0:j1) .* w(1:(j1-j0+1));
    end
end

% Raw RF plot
figure; 
plot(time * 1000000, signal2, 'b');
grid on; 
xlabel('Time (μs)',"FontSize", 16);
ylabel('Signal (V)',"FontSize", 16);
title('Raw RF Data over Time',"FontSize", 28);

%Verify FFt works
xw = signal2 .* hann(numel(signal2));
N  = numel(xw); df = fs/N;
f  = (-floor(N/2):ceil(N/2)-1).' * df;
X  = fftshift(fft(xw))/N;
magdB = 20*log10(abs(X)/max(abs(X)) + eps);
figure; plot(f/1e6, magdB); grid on; xlim([-30 30]); ylim([-100 0]);
xlabel('Frequency (MHz)'); ylabel('Mag (dB re. max)');
title('Windowed FFT of RF (post-blank)');

%% Gaussian Band Pass Filtering of Signal
t = time(:);
x = signal2(:);

% Sampling rate (assumes nearly-uniform sampling)
fs = 1 / mean(diff(t));

% Edit Based On Probe
f0 = 10e6; % center frequency (10 MHz probe)
FBW6  = 0.4; % fractional -6 dB bandwidth 
B6 = FBW6 * f0; % total -6 dB bandwidth

% Convert -6 dB half-width to Gaussian sigma:
sigma = (B6/2) / 1.177410022;

% Sanity Check (warn in Command Window)
f_low  = f0 - B6/2;
f_high = f0 + B6/2;
if fs < 2*f_high
    warning('Sampling rate %.2f MS/s may be too low for %.2f MHz upper band edge.', fs/1e6, f_high/1e6);
end

% Remove DC 
x = x - mean(x);

% FFT-domain zero-phase Gaussian BPF
N  = numel(x);
df = fs / N;
f = (-floor(N/2):ceil(N/2)-1).' * df;  % centered frequency grid (N points)
X = fftshift(fft(x));

% Plot the FFT
figure;
plot(f/1e6, abs(X), 'LineWidth', 2, 'Color', 'b');
xlabel('Frequency (MHz)', "FontSize", 16);
ylabel('Magnitude', "FontSize", 16);
title('FFT of the Raw Data', "FontSize", 28);
axis([-2500 2500 0 100])
grid on;

% Two Gaussian lobes at ±f0
Hpos = exp(-0.5*((f - f0)./sigma).^2);
Hneg = exp(-0.5*((f + f0)./sigma).^2);
H = Hpos + Hneg;
H = H / max(H); % unity gain in-band

% Plot the mask function H
figure;
plot(f/1e6, H, 'LineWidth', 3,'Color','b'); 
grid on; 
hold on;
xlabel('Frequency (MHz)',"FontSize",16);
ylabel('|H(f)|',"FontSize",16);
title(sprintf('Gaussian Bandpass (f0 = %.1f MHz, Bandwidth = %.0f%%)', f0/1e6, FBW6*100),"FontSize",16);
axis([0 20 0 1])

Y = X .* H;
%Plot the filter data in frequnecy domain
figure
plot(f/1e6, abs(Y), 'LineWidth', 2, 'Color', 'b');
xlabel('Frequency (MHz)', "FontSize", 16);
ylabel('Magnitude', "FontSize", 16);
title('Filtered FFT of the Raw Data', "FontSize", 28);
grid on;
axis([-50 50 0 150])

%Inverse Fourier Transform back to time domain
y = real(ifft(ifftshift(Y))); % filtered RF 

%% Analytic signal and envelope extraction with Hilbert Transform
xa  = hilbert(y);
envelopeSignal = abs(xa);

% Mild envelope smoothing (e.g., 0.5 us window)
envelopeSignal = movmean(envelopeSignal, 350);

% Convert time to depth
depth = t * 1540 / 2; % Assuming speed of sound in tissue is 1540 m/s
figure
plot(depth, y,'Color','b'); 
grid on;
xlabel('Depth (m)',"FontSize",16);
ylabel('Signal Amplitude (V)',"FontSize",16);
title(sprintf('Gaussian Bandpass Filtered RF Signal (f0=%.1f MHz, FBW6=%.0f%%)', f0/1e6, FBW6*100),'FontSize',28);

%% Analytic signal and envelope extraction with Hilbert Transform Without Function
filteredSignal = y(:); % real RF after your Gaussian BPF
N = numel(filteredSignal);

% FFT (centered) and magnitude plot
X = fftshift(fft(filteredSignal));
figure; 
plot(f/1e6, abs(X), 'LineWidth', 1.5);
xlabel('Frequency (MHz)',"FontSize",16);
ylabel('|X(f)|',"FontSize",16); 
title('FFT of Filtered Signal',"FontSize",28); 
axis([-20 20 0 150])
grid on;

XH = X;                               
% -90° on positive freqs, +90° on negative freqs
XH(f > 0) = XH(f > 0) * (-1i); % multiply by -j
figure; 
plot(f/1e6, angle(XH), 'LineWidth', 1.5);
title('test')
XH(f < 0) = XH(f < 0) * (+1i); % multiply by +j
figure; 
plot(f/1e6, angle(XH), 'LineWidth', 1.5);
title('test')

% DC
XH(f == 0) = 0;

%Phase Shift Transform
xhat = real(ifft(ifftshift(XH)));

% Analytic signal and envelope
xa = filteredSignal + 1i*xhat; % analytic signal
env  = abs(xa); % envelope 

% Plot RF and envelope together
figure; 
plot(depth, env, 'r', 'LineWidth', 1.5);
xlabel('Depth (m)'); 
ylabel('Amplitude'); 
grid on;
title('RF and Hilbert Envelope - No Function');

%% Looking at the 
X = fftshift(fft(filteredSignal));
XH = X; % apply Hilbert multiplier
XH(f>0) = XH(f>0) * (-1i); % −90° phase shift on +freqs
XH(f<0) = XH(f<0) * (+1i); % +90° phase shift on −freqs
XH(f==0) = 0;

phi0 = angle(X);
phi1 = angle(XH);
dphi = angle(exp(1i*(phi1 - phi0))); % wrapped Δphase in (−π,π]

mag  = abs(X);
mask = mag > 0.01*max(mag); % ignore tiny bins

figure; tiledlayout(3,1,'TileSpacing','compact');
nexttile; 
plot(f/1e6, mag); 
grid on;
xlabel('Frequency (MHz)'); 
ylabel('|X|'); 
title('Spectrum magnitude');

nexttile; 
plot(f(mask)/1e6, unwrap(phi0(mask)), 'b'); 
hold on;
plot(f(mask)/1e6, unwrap(phi1(mask)), 'r'); 
grid on;
xlabel('Frequency (MHz)'); 
ylabel('Phase (rad)');
legend('Before','After'); 
title('Phase before/after Hilbert');

nexttile; 
plot(f(mask)/1e6, dphi(mask)*180/pi, 'k'); 
hold on;
yline(-90,'--r','-90° (f>0)'); 
yline(+90,'--b','+90° (f<0)'); yline(0,'k:');
grid on; xlabel('Frequency (MHz)'); 
ylabel('\Delta phase (deg)');
title('Hilbert Δphase: −90° on +f, +90° on −f');


%% 
figure
plot(t, y,'Color','b'); 
grid on;
xlabel('Depth (m)',"FontSize",16);
ylabel('Signal Amplitude (V)',"FontSize",16);
title(sprintf('Gaussian Band Filtered RF Signal (f0=%.1f MHz, FBW6=%.0f%%)', f0/1e6, FBW6*100),'FontSize',28);

figure
plot(depth, envelopeSignal, 'LineWidth', 2,'Color','red');
ylabel("Signal Amplitude (V)","FontSize",16);
xlabel("Depth (m)","FontSize",16);
title('Envelope Signal','FontSize',28);
grid on;

figure
plot(depth, envelopeSignal, 'LineWidth', 2,'Color','red');
hold on
plot(depth,y,'Color','b')
ylabel("Signal Amplitude (V)","FontSize",16);
xlabel("Depth (m)","FontSize",16);
title('Envelope and Raw Signal','FontSize',28);
legend('Envelope','Raw')
grid on;


%% Log Compression and Normalization - Be Careful Doing This. Very Sample Dependent
e = max(envelopeSignal, 0); % guard tiny negatives
eps0 = 1e-12; % absolute floor (linear units)
pFloor = 3; % 1–5% typical
floorVal = max(eps0, prctile(e(:), pFloor));
e_dB = 20*log10(max(e, floorVal));

% Dynamic-range window(clip) 
DRdB = 50; %change
hi = max(e_dB(:));
lo = hi - DRdB;
logCompressedSignal = min(max(e_dB, lo), hi);

figure; 
plot(depth, logCompressedSignal,'Color','b'); 
ylabel("Log Voltage (dB)","FontSize",16);
xlabel("Depth (m)","FontSize",16);
title('Log Compressed Signal','FontSize',28); 
grid on;

% Normalize the log-compressed signal for better visualization
normalizedSignal = (logCompressedSignal - min(logCompressedSignal)) / (max(logCompressedSignal) - min(logCompressedSignal));
normalizedSignal(normalizedSignal < 0) = 0;

% Plot the normalized signal
figure;
plot(depth, normalizedSignal,'Color','b');
ylabel("Normalized Signal","FontSize",16);
xlabel("Depth (m)","FontSize",16);
title('Normalized Log Compressed Signal','FontSize',28);
grid on;


%% Ultrasound A-scan image (single column)
% Display as an image: one column  
gamma = 0.7;
normalizedSignal = normalizedSignal .^ gamma; %gamma scaling
z = normalizedSignal(:);
D = depth(:);                 

% Replicate to 2 columns so it's visible
C = [z z];

figure
imagesc(1, D*1000, C); % depth in mm;
colorbar;
colormap('gray')
xlabel('A-line index',"FontSize",16);
ylabel('Depth (mm)',"FontSize",16);
title('Ultrasound A-scan as Image Column','FontSize',28);
set(gca,'YDir','reverse'); 
