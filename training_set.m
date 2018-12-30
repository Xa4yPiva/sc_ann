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
signals = [xAM; xDSB; xLSB; xUSB; xFM];
sigsNum = size(signals, 1);
% modNames = ["AM", "DSB", "LSB", "USB", "FM", "Noise"];
% modNames = ["AM", "DSB", "FM", "Noise"];
modNames = ["AM", "DSB", "LSB", "USB", "FM"];
modNums = 1 : sigsNum;
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
lenEnv = size(envelopes, 2);
N = 2 ^ nextpow2(lenEnv);
freqs = (-N/2 : (N-1)/2) * fs2 / N;
specEnv = fftshift(abs(fft(envelopes, N, 2)) / (N/2), 2);

%% Plot signals
% figure(1);
% subplot(5,1,1); plot(freqs, specAM); grid on;
% subplot(5,1,2); plot(freqs, specDSB); grid on;
% subplot(5,1,3); plot(freqs, specLSB); grid on;
% subplot(5,1,4); plot(freqs, specUSB); grid on;
% subplot(5,1,5); plot(freqs, specFM); grid on;

% figure(2);
% subplot(5,1,1); plot(freqs, mag2db(specEnv(1,:))); grid on;
% subplot(5,1,2); plot(freqs, mag2db(specEnv(2,:))); grid on;
% subplot(5,1,3); plot(freqs, mag2db(specEnv(3,:))); grid on;
% subplot(5,1,4); plot(freqs, mag2db(specEnv(4,:))); grid on;
% subplot(5,1,5); plot(freqs, mag2db(specEnv(5,:))); grid on;
% return
%% KFs
lenPowMin = 12;
lenPowMax = floor(log2(lenSignal));
lenFrameMin = 2 ^ lenPowMin;
lenFrameMax = 2 ^ lenPowMax;
thresholds.ampl = 1;
expNum = 10000;
% lenFrame = zeros(1, sigsNum * expNum);
lenFrame = 2^12;
idxF = 1;
cyclesNum = sigsNum * expNum;
% h = waitbar(0, 'Computing...');
disp("Computing training set ...");
tic
for j = 1 : expNum
    snr = rand() * 30 - 15;
    pos = randi([0, lenEnv - lenFrame-1]);
    for i = 1 : sigsNum
%         lenFrame(idxF) = 2 ^ (round(rand() * (lenPowMax - lenPowMin) + lenPowMin));
%         pos = round(rand() * (lenSignal - lenFrame));
        env = awgn(envelopes(i, pos+1 : pos+lenFrame), snr, 'measured');
        kf = KeyFeatures(env, thresholds.ampl);
%         if (~isreal([kf.gammaMax, kf.sigmaDP, kf.sigmaAP]))
%             kf = KeyFeatures(env, thresholds.ampl);
%         end
        ts(idxF).kf = kf;
%         ts(idxF).lenFrame = lenFrame(idxF);
        ts(idxF).lenFrame = lenFrame;
        ts(idxF).snr = snr;
%         ts(idxF).type = modNames(i);
        ts(idxF).type = modNums(i);
%         waitbar(idxF / cyclesNum);
        idxF = idxF + 1;
    end
end
toc
disp('... done.')
% close(h);

kfs = [ts(:).kf];
ts_inputs = [abs([[kfs.gammaMax]; [kfs.sigmaDP]; [kfs.sigmaAP]]); [kfs.P]];
ts_targets = repmat(diag(ones(1, sigsNum)), 1, expNum);











