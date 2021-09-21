clc
clear
close all

load example.mat % 3x 4030 signals of EEG, EOG1, EOG2
% We will only use EEG signal for now
% t = 0:0.01:10;
fs = 1000;
fn = fs/2;
eeg = example(1,:);
% eeg = 10*sin(2*pi*0.5*t) + 8*sin(2*pi*5*t) + 7*sin(2*pi*7*t) + 4*sin(2*pi*12*t) + 3.5*sin(2*pi*15*t) + 2*sin(2*pi*20*t) + 1.5*sin(2*pi*50*t) + sin(2*pi*100*t);
offset = mean(eeg);
eeg = eeg - offset;
%Check your signal
figure
plot(eeg)
t = 0:1/fs:(length(eeg)-1)/fs;
%FFT

Y = fft(eeg);
L = numel(Y);
P2 = abs(Y / L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
df = fs/length(eeg);
f = 0:df:fs/2;

% plot signal in frequency domain
figure
plot(f, P1)
xlim([0,100])
% seems like there is a strong motion artifact/noise

% Filter noise
% fc = 0.01;
% filtered_eeg = highpass(eeg, fc, fs);

figure 
plot(eeg)


% get EEG bandwidth for delta, theta, alpha, beta
filter_order = 3;
% delta wave (1-4Hz)
[z, p, k] = butter(filter_order, [1/fn 4/fn], "bandpass");
[b,a] = zp2tf(z,p,k);
delta = filtfilt(b,a,eeg);

%plot delta
figure
plot(t,delta)
title("Delta (1-4Hz)")


%alpha (8-12)
[z,p,k] = butter(filter_order, [8/fn 12/fn], "bandpass");
[b,a] = zp2tf(z,p,k);
alpha = filtfilt(b,a,eeg);

%plot alpha
figure
plot(t,alpha)
title("Alpha(8-12Hz)")

%theta (4-8)
[z, p, k] = butter(filter_order, [4/fn 8/fn], "bandpass");
[b,a] = zp2tf(z,p,k);
theta = filtfilt(b,a,eeg);

%plot theta
figure
plot(t,theta)
title("Theta(4-8Hz)")
% 
%beta (12-35)
[z, p, k] = butter(filter_order, [12/fn, 35/fn], "bandpass");
[b,a] = zp2tf(z,p,k);
beta = filtfilt(b,a,eeg);

%plot beta
figure
plot(t,beta)
title("Beta(12-35Hz)")



