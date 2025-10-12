close all
clear all

%Import CSV File
rawRFData = readtable("tek0001.csv");
clc

%% Extract relevant columns from the imported data
time = rawRFData.WaveformType; % Assuming 'Time' is a column in the CSV
signal = rawRFData.ANALOG; % Assuming 'Signal' is another column
synch = rawRFData.Var3;

%Remove the transient spike from the probe
time(1:400000) = []; 
signal(1:400000) = []; 
synch(1:400000) = [];

%Set time to start at 0
time(:) = time(:) - min(time(:));

%Raw RF over time
plot(time, signal, 'Color','blue');
xlabel('Time (s)');
ylabel('Signal Amplitude');
title('Raw RF Signal Over Time');
grid on;

%% Gaussian Band Pass Filtering of Signal
t = time(:);
x = signal(:);

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

% Remove DC (helps analytic signal)
x = x - mean(x);

% FFT-domain zero-phase Gaussian BPF
N  = numel(x);
df = fs / N;
f  = (-floor(N/2):ceil(N/2)-1).' * df;  % centered frequency grid (N points)
X = fftshift(fft(x));

% Two Gaussian lobes at ±f0
Hpos = exp(-0.5*((f - f0)./sigma).^2);
Hneg = exp(-0.5*((f + f0)./sigma).^2);
H = Hpos + Hneg;
H = H / max(H); % unity gain in-band

Y = X .* H;
y = real(ifft(ifftshift(Y))); % filtered RF (zero-phase)

% Analytic signal and envelope (pre-Hilbert filter already applied)
xa  = hilbert(y);
envelopeSignal = abs(xa);

% Mild envelope smoothing (e.g., 0.5 us window)
envelopeSignal = movmean(envelopeSignal, 350);

% Convert time to depth
depth = t * 1540 / 2; % Assuming speed of sound in tissue is 1540 m/s

% Plots
figure
plot(depth, y); 
grid on;
xlabel('Depth (m)'); 
ylabel('Amplitude'); 
title(sprintf('Gaussian Band Filtered RF Signal (f0=%.1f MHz, FBW6=%.0f%%)', f0/1e6, FBW6*100));

figure
plot(depth, envelopeSignal, 'LineWidth', 2,'Color','red');
ylabel("Envelope (V)"); 
xlabel("Depth (m)");
title('Envelope Signal');
grid on;

figure
plot(depth, envelopeSignal, 'LineWidth', 2,'Color','red');
hold on
plot(depth,y,'Color','b')
ylabel("Envelope (V)"); 
xlabel("Depth (m)");
title('Envelope and Raw Signal');
legend('Envelope','Raw')
grid on;


%% Log Compression and Normalization
% Floor to avoid log(0) and set background level
e = max(envelopeSignal, 0); % guard tiny negatives
eps0 = 1e-12; % absolute floor (linear units)
pFloor = 3; % 1–5% typical
floorVal = max(eps0, prctile(e(:), pFloor));

e_dB = 20*log10(max(e, floorVal));

% Dynamic-range window(clip) 
DRdB = 45; %change
hi = max(e_dB(:));
lo = hi - DRdB;
logCompressedSignal = min(max(e_dB, lo), hi);

figure; 
plot(depth, logCompressedSignal); 
ylabel("Log Voltage (dB)"); 
xlabel("Depth (m)"); 
title('Log Compressed Signal'); 
grid on;

% Normalize the log-compressed signal for better visualization
normalizedSignal = (logCompressedSignal - min(logCompressedSignal)) / (max(logCompressedSignal) - min(logCompressedSignal));
normalizedSignal(normalizedSignal < 0) = 0;

% Plot the normalized signal
figure;
plot(depth, normalizedSignal);
ylabel("Normalized Signal");
xlabel("Depth (m)");
title('Normalized Log Compressed Signal');
grid on;


%% Ultrasound A-scan image (single column)
% Display as an image: one column
gamma = 0.8;                         
normalizedSignal = normalizedSignal .^ gamma;
z = normalizedSignal(:);
D = depth(:);                 

% Replicate to 2 columns so it's visible
C = [z z];

figure
imagesc([1 2], D*1000, C);% depth in mm; x=two columns  
colorbar;
colormap('gray')
xlabel('A-line index');
ylabel('Depth (mm)');
title('Ultrasound A-scan as Image Column');
set(gca,'YDir','reverse'); 
