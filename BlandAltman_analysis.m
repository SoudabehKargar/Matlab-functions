clc
clear variables
close all


load('Exam_5633_5645_ROIvalues_all.mat')

%% Ktrans
close all
fig = figure;

subplot(1,4,1)
data1 = roiMean_Ktrans_1;
data2 = roiMean_Ktrans_2;
[means,diffs,meanDiff,CR,linFit] = BlandAltman(data1,data2,2)
text(min(means)*1.05,CR(1)*.9,num2str(CR(1)),'fontsize',14,'color','b')
text(min(means)*1.05,CR(2)*.75,num2str(CR(2)),'fontsize',14,'color','b')
text(min(means)*1.05,CR(1)*1.15,'+1.96\sigma','fontsize',14,'color','r')
text(min(means)*1.05,CR(2)*1.25,'-1.96\sigma','fontsize',14,'color','r')
ylim([-.2 .3])
title('K^{trans} [min^{-1}]')

subplot(1,4,2)
data1 = roiMean_vp_1;
data2 = roiMean_vp_2;
[means,diffs,meanDiff,CR,linFit] = BlandAltman(data1,data2,2)
text(min(means)*1.05,CR(1)*.9,num2str(CR(1)),'fontsize',14,'color','b')
text(min(means)*1.05,CR(2)*.8,num2str(CR(2)),'fontsize',14,'color','b')
text(min(means)*1.05,CR(1)*1.1,'+1.96\sigma','fontsize',14,'color','r')
text(min(means)*1.05,CR(2)*1.15,'-1.96\sigma','fontsize',14,'color','r')
ylim([-.07 .09])
xlim([0.03 .1])
title('v_p [a. u.]')

subplot(1,4,3)
data1 = roiMean_ve_1;
data2 = roiMean_ve_2;
[means,diffs,meanDiff,CR,linFit] = BlandAltman(data1,data2,2)
text(min(means)*1.05,CR(1)*.8,num2str(CR(1)),'fontsize',14,'color','b')
text(min(means)*1.05,CR(2)*.8,num2str(CR(2)),'fontsize',14,'color','b')
text(min(means)*1.05,CR(1)*1.2,'+1.96\sigma','fontsize',14,'color','r')
text(min(means)*1.05,CR(2)*1.2,'-1.96\sigma','fontsize',14,'color','r')
ylim([-.18 .1])
title('v_e [a. u.]')

subplot(1,4,4)
data1 = roiMean_T10_1;
data2 = roiMean_T10_2;
[means,diffs,meanDiff,CR,linFit] = BlandAltman(data1,data2,2)
text(min(means)*1.05,CR(1)*.85,num2str(CR(1)),'fontsize',14,'color','b')
text(min(means)*1.05,CR(2)*.85,num2str(CR(2)),'fontsize',14,'color','b')
text(min(means)*1.05,CR(1)*1.15,'+1.96\sigma','fontsize',14,'color','r')
text(min(means)*1.05,CR(2)*1.15,'-1.96\sigma','fontsize',14,'color','r')
ylim([-350 270])
title('T_{1} [msec]')
set(fig,'Units', 'normalized', 'Position', [0.2,0.2,.6,.2])

