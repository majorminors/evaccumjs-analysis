function betterThresholdPlotting(res1, res2)

plotOptions.CIthresh       = true;

plotOptions2 = plotOptions;
plotOptions2.lineColor = [17,17,17]/255; % plot grey because it should be the same, but if not de-emphasise it

[hline,hdata] = plotPsych(res1,plotOptions);

hold on

[hline2,hdata2] = plotPsych(res2,plotOptions2);

%legend([hline,hline2],'Fit1','Fit2')

return
end