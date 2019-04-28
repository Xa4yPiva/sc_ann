function [] = PlotManyProbs(p1, p2, snr1, snr2, markers, names)

figsNum = size(p, 3);
colors = lines(max(min(size(p1), min(size(p2)))));
figure(7);
for i = 1 : sigsNum
%     plot(snrANN, pANN(i,:), 'marker', mType(i), 'markersize', 10, 'linewidth', 2, 'color', colors(i,:));
    plot(snr1, p1(i,:), 'marker', markers(i), 'linewidth', 2, 'color', colors(i,:));
    hold on;
end
for i = 1 : sigsNum
    m = find(ld.decRight == modNames(i));
%     plot(snrT, pT(m,:), '--', 'marker', mType(i), 'markersize', 10, 'linewidth', 2, 'color', colors(i,:));
    plot(snr2, p2(m,:), '--', 'marker', markers(i), 'linewidth', 2, 'color', colors(i,:));
    hold on;
end
mods = [strcat(modNames, "-ANN"), strcat(modNames, "-Thresholds")];
grid on;
title('Thresholds vs ANN'); xlabel('SNR, dB'); ylabel('Probability of right decision');
legend(mods, 'location', 'southeast'); legend('show');
set(gcf, 'color', 'w'); set(groot, 'DefaultAxesFontSize', 18);

end

