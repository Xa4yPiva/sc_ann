function [ pRight pRightOverall ] = ProbRightDecision( net, envelopes, thresholds, decisionsRight, SNR, expNum, lenFrame )
%PROBRIGHTDECISION Summary of this function goes here
%   Detailed explanation goes here

sigsNum = min(size(envelopes));
framesNum = floor(max(size(envelopes)) / lenFrame);
lenSNR = length(SNR);
pRight = zeros(sigsNum, lenSNR);
pRightOverall = zeros(1, lenSNR);
lenEnv = size(envelopes, 2);
% cyclesNum = sigsNum * lenSNR;
% iteration = 0;
% h = waitbar(0, 'Computing probability of right decision...');
disp('Computing probability of right decision ...');
tic
parfor i = 1 : lenSNR
    decRightOverallNum = 0;
    for k = 1 : sigsNum
        decRightNum = 0;
        for j = 1 : expNum
            pos = randi([1, lenEnv - lenFrame]);
            env = awgn(envelopes(k, pos : pos+lenFrame), SNR(i), 'measured');
%             env = awgn(envelopes(k, pos*lenFrame+1 : (pos+1)*lenFrame), SNR(i), 'measured');
            kf = KeyFeatures(env, thresholds.ampl);
            netInputs = [abs([kf.gammaMax; kf.sigmaDP; kf.sigmaAP]); kf.P];
            netOutputs = net(netInputs);
            [oMax, oIdx] = max(netOutputs);
%             decision = decisionsRight(oIdx);
            if (oIdx == k)
                decRightNum = decRightNum + 1;
            end
        end
        decRightOverallNum = decRightOverallNum + decRightNum;
        pRight(k, i) = decRightNum / expNum;
%         iteration = iteration + 1;
%         waitbar(iteration / cyclesNum);
    end
    pRightOverall(i) = decRightOverallNum / (sigsNum * expNum);
end
toc
% close(h);
disp('... done.');

end
