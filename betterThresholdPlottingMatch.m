function betterThresholdPlottingMatch(res1, res2)
% will plot tmp1 green and tmp2 red

plotOptions.CIthresh       = true;


plotOptions1 = plotOptions;
plotOptions1.lineColor = [0,1,0];
plotOptions1.dataColor      = [0,1,0];

plotOptions2 = plotOptions;
plotOptions2.lineColor = [1,0,0];
plotOptions2.dataColor      = [1,0,0];

[hline,hdata] = plotPsych(res1,plotOptions1);

hold on

[hline2,hdata2] = plotPsych(res2,plotOptions2);

legend([hline,hline2],'easy dots','hard dots')

return
end