addpath(genpath('../matlab_utils'));
% addpath(genpath('../sc_amra'));
addpath(genpath('../sc_am_fm'));
addpath(genpath('../sc_common'));
clear;
close all;

%% Signal
fs0 = 12e3;
T = 1;
fLowHz = 300;
fHighHz = 3600;
bandHz = fHighHz - fLowHz;
data = RandomBandLimitedSignal(fs0, T, 20, fLowHz, fHighHz, 4000, 60, 1, 60, 'uniform');
fs = 200e3;
factor = fs / fs0;
[p, q] = rat(factor);
x = resample(data, p, q);
lenSignal = length(x);
N = 2 ^ nextpow2(lenSignal);
freqs = (-N/2 : (N-1)/2) * fs / N;
fc = 60e3;
offset = exp(1i * 2*pi*fc * (0:lenSignal-1)/fs);

%% Modulation
mAM = 0.3;
xAM = ammod(mAM * x, fc, fs, 0, 0.5);
specAM = fftshift(abs(fft(xAM, N))) / (N/2);
xDSB = ammod(x, fc, fs, 0, 0);
specDSB = fftshift(abs(fft(xDSB, N))) / (N/2);
xLSB = ssbmod(x, fc, fs, 0);
specLSB = fftshift(abs(fft(xLSB, N))) / (N/2);
xUSB = ssbmod(x, fc, fs, 0, 'upper');
specUSB = fftshift(abs(fft(xUSB, N))) / (N/2);
fDev = 25e3;
xFM = fmmod(x, fc, fs, fDev);
specFM = fftshift(abs(fft(xFM, N))) / (N/2);

xNoise = zeros(1, lenSignal);


%% Signals, Decisions
% signals = [xAM; xDSB; xLSB; xUSB; xFM; xNoise];
% signals = [xAM; xDSB; xFM; xNoise];
% signals = [xAM; xDSB; xFM;];
signals = [xAM; xDSB; xLSB; xUSB; xFM;];
sigsNum = size(signals, 1);
% modNames = ["AM", "DSB", "LSB", "USB", "FM", "Noise"];
% modNames = ["AM", "DSB", "FM", "Noise"];
% modNames = ["AM", "DSB", "FM"];
modNames = ["AM", "DSB", "LSB", "USB", "FM"];
modNums = 1 : sigsNum;
mType = ['x', '*', 'o', 'd', 's', '+'];

envelopes = zeros(size(signals));
% offsets = [-fc, -fc, -fc + bandHz/2 + fLowHz, -fc - bandHz/2 - fLowHz, -fc, 0];
offsets = -fc * ones(1, sigsNum);
expOff = exp(1i * 2*pi*offsets' .* (0:lenSignal-1)/fs);
for i = 1 : sigsNum
    envelopes(i, :) = hilbert(signals(i, :)) .* expOff(i, :);
end
fs2 = 80e3;
factor = fs2 / fs;
[p, q] = rat(factor);
envelopes = (resample(envelopes', p, q))';

% specEnv = fftshift(abs(fft(envelopes, N, 2)) / (N/2), 2);

%% Probabilities of Right Decision
% load('network1.mat');
% load('network2_lf8k.mat');
% load('network4_3-7-10-4_logsig_lf8k.mat');
% load('network5_3-7-10-7-4_tansig_lin_lf8k.mat');
% load('network6_3-7-10-7-6_tansig_lf8k.mat');
% load('network7_3-7-10-7-3_tansig_lf8k.mat');
load('net4-7-5.mat');
net = net4_7_5;
thresholds.ampl = 1;
snr = -15 : 2 : 10;
expNum = 1e3;
lenFrame = 2^12;
% [pRight, pro] = ProbRightDecision(net, envelopes, thresholds, modNames, snr, expNum, lenFrame);

%% Plot pRight
close all;
tNames = [];
PlotProbRight([6,7], pRight, pro, snr, mType, modNames, tNames);

%% Plot pErr
pErr = 1 - pRight;
peo = 1 - pro;
PlotProbError([8,9], pErr, peo, snr, mType, modNames, tNames);

%% Plot Probabilities Thresholds vs ANN
close all;
pANN = pRight;
snrANN = snr;
ld = load('../sc_amra/am_dsb_lsb_usb_fm.mat', 'pRight', 'snr', 'decRight');
pT = ld.pRight;
snrT = ld.snr;
colors = lines(sigsNum);
figure(7);
for i = 1 : sigsNum
%     plot(snrANN, pANN(i,:), 'marker', mType(i), 'markersize', 10, 'linewidth', 2, 'color', colors(i,:));
    plot(snrANN, pANN(i,:), 'marker', mType(i), 'linewidth', 2, 'color', colors(i,:));
    hold on;
end
for i = 1 : sigsNum
    m = find(ld.decRight == modNames(i));
%     plot(snrT, pT(m,:), '--', 'marker', mType(i), 'markersize', 10, 'linewidth', 2, 'color', colors(i,:));
    plot(snrT, pT(m,:), '--', 'marker', mType(i), 'linewidth', 2, 'color', colors(i,:));
    hold on;
end
mods = [strcat(modNames, "-ANN"), strcat(modNames, "-Thresholds")];
grid on;
title('Thresholds vs ANN'); xlabel('SNR, dB'); ylabel('Probability of right decision');
legend(mods, 'location', 'southeast'); legend('show');
set(gcf, 'color', 'w'); set(groot, 'DefaultAxesFontSize', 18);

%% Plot Probabilities of error T vs ANN
close all;
pErrANN = 1 - pRight;
snrANN = snr;
ld = load('../sc_amra/am_dsb_lsb_usb_fm.mat', 'pRight', 'snr', 'decRight');
pErrT = 1 - ld.pRight;
snrT = ld.snr;
colors = lines(sigsNum);
figure(8);
for i = 1 : sigsNum
    semilogy(snrANN, pErrANN(i,:), 'marker', mType(i), 'linewidth', 2, 'color', colors(i,:));
    hold on;
end
for i = 1 : sigsNum
    m = find(ld.decRight == modNames(i));
    semilogy(snrT, pErrT(m,:), '--', 'marker', mType(i), 'linewidth', 2, 'color', colors(i,:));
    hold on;
end
mods = [strcat(modNames, "-ANN"), strcat(modNames, "-Thresholds")];
grid on;
title('Thresholds vs ANN'); xlabel('SNR, dB'); ylabel('Probability of error');
legend(mods, 'location', 'northeast'); legend('show');
set(gcf, 'color', 'w'); set(groot, 'DefaultAxesFontSize', 18);
