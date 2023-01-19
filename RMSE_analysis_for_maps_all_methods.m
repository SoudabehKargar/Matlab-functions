clc
clear variables
close all

%% 
filter_low_signal = 'n';
cutoff = .2;
if filter_low_signal == 'n'
    cutoff = '';
end

SNR = [1000;200;100;90;80;70;60;50;40;30;20;10;5];
KtransCaxis = [0,1.2];
vpCaxis     = [0,.3];
veCaxis     = [0,1.4];
kepCaxis    = [0,2];
T10Caxis    = [0,4000];
S0Caxis     = [0 0.6];
BATCaxis    = [0 40];
dictionary_size = '200k'
for ii=1:3
    if dictionary_size == '400k'
        if ii==1
            % StepLike  282/280  289/290
            version = 55            % 48/51/55
            model_1 = 290
            model_2 = 289
            AIF_type = 'StepLike'

        elseif ii==2
            % Parker    254/256 291/292
            version = 53        % 47/50/53
            model_1 = 292
            model_2 = 291
            AIF_type = 'OriginalParker'

        elseif ii==3
            % Modified Parker   275/281 287/288
            version = 54        % 49/52/54
            model_1 = 288
            model_2 = 287
            AIF_type = 'ModifiedParker'
        end
    elseif dictionary_size == '200k'
        if ii==1
            % StepLike  282/280  289/290
            version = 56        % 48/51/55
            model_1 = 296
            model_2 = 295
            AIF_type = 'StepLike'

        elseif ii==2
            % Parker    254/256 291/292
            version = 53        % 47/50/53
            model_1 = 298
            model_2 = 297
            AIF_type = 'OriginalParker'

        elseif ii==3
            % Modified Parker   275/281 287/288
            version = 54        % 49/52/54
            model_1 = 294
            model_2 = 293
            AIF_type = 'ModifiedParker'
        end
    end
    
    
filename5 = ['NRMSE_analysis_Version_',num2str(version),'_AIFType_',num2str(AIF_type),'_model_',num2str(model_1),'_from_',num2str(model_2),'_Cutoff_',num2str(filter_low_signal),'_',num2str(cutoff),'.mat']
recon5 = load(filename5);

% save(filename5, 'NRMSE_Ktrans','RMSE_Ktrans','NRMSE_vp','RMSE_vp','NRMSE_ve','RMSE_ve','NRMSE_T10','RMSE_T10','NRMSE_S0','RMSE_S0','NRMSE_BAT','RMSE_BAT','SNR','RMSE_fit_mean')

NRMSE_Ktrans_all(ii,:)  = recon5.NRMSE_Ktrans;
NRMSE_vp_all(ii,:)      = recon5.NRMSE_vp;
NRMSE_ve_all(ii,:)      = recon5.NRMSE_ve;
NRMSE_T10_all(ii,:)     = recon5.NRMSE_T10;
NRMSE_S0_all(ii,:)      = recon5.NRMSE_S0;
NRMSE_BAT_all(ii,:)     = recon5.NRMSE_BAT;

RMSE_fit_mean_all(ii,:) = recon5.RMSE_fit_mean;
str{ii} = AIF_type; 
str2{ii} = version;
str3{ii} = model_1;
str4{ii} = model_2;
str5{ii} = [str{ii},': ','v.',num2str(str2{ii}),' model:',num2str(str3{ii}),'/',num2str(str4{ii})]
end
%%
% close all
fig = figure
x = [100 ;90; 80; 70 ;60 ;50 ;40; 30; 20 ]
y = [RMSE_fit_mean_all(1,3:end-2)', RMSE_fit_mean_all(2,3:end-2)', RMSE_fit_mean_all(3,3:end-2)'];
b = bar(x,y)
legend(str)
text(70,.8*max(y(:)),str5,'Color','k','FontSize',14)
xlabel('SNR')
ylabel('Mean Square Error (Goodness of fit)')
title(['Dictionary Size: ',num2str(dictionary_size),' - Cutoff:',num2str(cutoff)])
set(gca,'FontSize',14)
set(fig,'Units', 'normalized', 'Position', [0.2,0.3,.45,.45])


fig = figure
x = [100 ;90; 80; 70 ;60 ;50 ;40; 30; 20 ;10; ]
y = [RMSE_fit_mean_all(1,3:end-1)', RMSE_fit_mean_all(2,3:end-1)', RMSE_fit_mean_all(3,3:end-1)'];
b = bar(x,y)
legend(str)
text(70,.8*max(y(:)),str5,'Color','k','FontSize',14)
xlabel('SNR')
% ylim([0 3.5])
title(['Dictionary Size: ',num2str(dictionary_size),' - Cutoff:',num2str(cutoff)])
ylabel('Mean Square Error (Goodness of fit)')
set(gca,'FontSize',14)
set(fig,'Units', 'normalized', 'Position', [0.2,0.3,.45,.45])

fig = figure
x = [200; 100 ;90; 80; 70 ;60 ;50 ;40; 30; 20 ;10; ]
y = [RMSE_fit_mean_all(1,2:end-1)', RMSE_fit_mean_all(2,2:end-1)', RMSE_fit_mean_all(3,2:end-1)'];
b = bar(x,y)
legend(str)
text(70,.8*max(y(:)),str5,'Color','k','FontSize',14)
xlabel('SNR')
title(['Dictionary Size: ',num2str(dictionary_size),' - Cutoff:',num2str(cutoff)])
ylabel('Mean Square Error (Goodness of fit)')
set(gca,'FontSize',14)
set(fig,'Units', 'normalized', 'Position', [0.2,0.3,.45,.45])

fig = figure
x = [1000; 200; 100 ;90; 80; 70 ;60 ;50 ;40; 30; 20 ;10; ]
y = [RMSE_fit_mean_all(1,1:end-1)', RMSE_fit_mean_all(2,1:end-1)', RMSE_fit_mean_all(3,1:end-1)'];
b = bar(x,y)
legend(str)
text(70,.8*max(y(:)),str5,'Color','k','FontSize',14)
xlabel('SNR')
title(['Dictionary Size: ',num2str(dictionary_size),' - Cutoff:',num2str(cutoff)])
ylabel('Mean Square Error (Goodness of fit)')
set(gca,'FontSize',14)
set(fig,'Units', 'normalized', 'Position', [0.2,0.3,.45,.45])

% xtips1 = b(1).XEndPoints
% ytips1 = b(1).YEndPoints
% labels1 = string(b(1).YData)
% text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
%     'VerticalAlignment','bottom')
% 
% xtips2 = b(2).XEndPoints;
% ytips2 = b(2).YEndPoints;
% labels2 = string(b(2).YData);
% text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
%     'VerticalAlignment','bottom')
% 
% xtips3 = b(3).XEndPoints;
% ytips3 = b(3).YEndPoints;
% labels3 = string(b(3).YData);
% text(xtips3,ytips3,labels3,'HorizontalAlignment','center',...
%     'VerticalAlignment','bottom')

%%
% StepLike  282/280  289/290
%         version = 51
%         model_1 = 290
%         model_2 = 289
%         AIF_type = 'StepLike'
%     
%         % Parker    254/256 291/292
%         version = 50
%         model_1 = 292
%         model_2 = 291
%         AIF_type = 'OriginalParker'
%     
%         % Modified Parker   275/281 287/288
%         version = 52
%         model_1 = 288
%         model_2 = 287
%         AIF_type = 'ModifiedParker'

linestyle = {'-','--','-.'}
linestyle1 = {'_____','- - - -','- . - . -'}

fig = figure;
for ii=1:3
    plot(SNR,NRMSE_Ktrans_all(ii,:),num2str(linestyle{ii}),'marker','o','color',[0 0.4470 0.7410],'markerfacecolor',[0 0.4470 0.7410],'linewidth',1.5), hold on 
    plot(SNR,NRMSE_vp_all(ii,:),num2str(linestyle{ii}),'marker','*','color',[0.8500 0.3250 0.0980],'markerfacecolor',[0.8500 0.3250 0.0980],'linewidth',1.5), hold on 
    plot(SNR,NRMSE_ve_all(ii,:),num2str(linestyle{ii}),'marker','<','color',[0.9290 0.6940 0.1250],'markerfacecolor',[0.9290 0.6940 0.1250]	,'linewidth',1.5), hold on 
    plot(SNR,NRMSE_T10_all(ii,:),num2str(linestyle{ii}),'marker','x','color',[0.4940 0.1840 0.5560],'markerfacecolor',[0.4940 0.1840 0.5560],'linewidth',1.5), hold on 
    plot(SNR,NRMSE_S0_all(ii,:),num2str(linestyle{ii}),'marker','d','color',[0.4660 0.6740 0.1880],'markerfacecolor',[0.4660 0.6740 0.1880],'linewidth',1.5), hold on 
    plot(SNR,NRMSE_BAT_all(ii,:),num2str(linestyle{ii}),'marker','^','color',[0.3010 0.7450 0.9330],'markerfacecolor',[0.3010 0.7450 0.9330],'linewidth',1.5), hold on
    legend('Ktrans','vp','ve','T10','S0','BAT','linewidth',1.5)
    str1{ii} = [linestyle1{ii},'  ',str{ii},': ','v.',num2str(str2{ii}),' model:',num2str(str3{ii}),'/',num2str(str4{ii})]
end

xlabel('SNR')
ylabel('nRMSE (%)')
text(50,100,str1,'Color','k','FontSize',14)
title(['nRMSE Analysis Numeric Sim - Cutoff: ',num2str(filter_low_signal),' ',num2str(cutoff),'Dictionary Size: ',num2str(dictionary_size)])
xlim([0 100])
set(gca,'FontSize',14)
set(fig,'Units', 'normalized', 'Position', [0.2,0.3,.45,.45])

%% 
% clc
% A = randi(20,1,10)
% [B,I] = sort(A)
% 
% C = A(I)


    





