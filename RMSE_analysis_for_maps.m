%% Analyze RMSE vs SNR
clc
clear variables
close all

filter_low_signal = 'n';
cutoff = .2;
if filter_low_signal == 'n'
    cutoff = '';
end
plot_all_figures = 'n'
dictionary_size = '200k'
AIF = 'P';             %  P: Parker M: Modified Parker S: StepLike

filter_low_T1   = 'n';
T1_treshold     = 1000;

filter_large_error = 'y';
param_cutoff = 'T10';
param_error_cutoff = 20;

if dictionary_size == '400k'    
    if AIF == 'S'      % 282/280  289/290
        version = 51
        model_1 = 290
        model_2 = 289
        AIF_type = 'StepLike'

    elseif AIF == 'P'  % 254/256 291/292
        version = 50   % 47/50/53
        model_1 = 292
        model_2 = 291
        AIF_type = 'OriginalParker'

    elseif AIF == 'M'  % 275/281 287/288
        version = 52   % 49/52/54
        model_1 = 288
        model_2 = 287
        AIF_type = 'ModifiedParker'
    end
elseif dictionary_size == '200k'
    if AIF == 'S'      % 282/280 295/296
        version = 56   % 48/51/55
        model_1 = 296
        model_2 = 295
        AIF_type = 'StepLike'

    elseif AIF == 'P'  % 254/256 297/298
        version = 53
        model_1 = 298
        model_2 = 297
        AIF_type = 'OriginalParker'

    elseif AIF == 'M'  %  275/281 293/294
        version = 54
        model_1 = 294
        model_2 = 293
        AIF_type = 'ModifiedParker'
    end
end

%% 
SNR = [1000;200;100;90;80;70;60;50;40;30;20;10;5];
SNR_extra_dict = [1000];
%%
ktransCaxis = [0,1];
vpCaxis     = [0,.1];
veCaxis     = [0,1];
kepCaxis    = [0,1];
T10Caxis    = [0,2500];
S0Caxis     = [0 .5];
BATCaxis    = [0 30];
fCaxis      = [0 .5];
B1Caxis     = [0 1];
a1Caxis     = [0 10];
Error_max   = [0 20];

KtransCaxis = [0,1];
vpCaxis     = [0,.15];
veCaxis     = [0,1];
kepCaxis    = [0,2];
T10Caxis    = [0,3000];
S0Caxis     = [0.2 0.6];
BATCaxis    = [0 30];

KtransErrorCaxis = [0 100];
vpErrorCaxis = [0 100];
veErrorCaxis = [0 150];
T10ErrorCaxis = [0 100];
BATErrorCaxis = [0 20];
S0ErrorCaxis = [0 20];

%%
for ii=1:length(SNR)
    ii
    filename1 = ['/Volumes/MRIClinical/kargar/DL/DCE_DRONE_code/test_data_recon/test_data_recon_num_estimate_test_vals_v',num2str(version),'_SNR',num2str(SNR(ii)),'_model_',num2str(model_2),'_from_.mat'];
    recon1 = load(filename1);
    filename2 = ['/Volumes/MRIClinical/kargar/DL/DCE_DRONE_code/test_data_recon/test_data_recon_num_estimate_test_vals_v',num2str(version),'_SNR',num2str(SNR(ii)),'_model_',num2str(model_1),'_from_',num2str(model_2),'.mat'];
    recon2 = load(filename2);
    filename3 = ['/Volumes/MRIClinical/kargar/DL/DCE_DRONE_code/test_data_vals/test_data_vals_v',num2str(version),'.mat'];
    recon3 = load(filename3);
    filename4 = ['/Volumes/MRIClinical/kargar/DL/DCE_DRONE_code/data_to_match/data_to_match_recon_estimate_num_test_vals_v',num2str(version),'_SNR',num2str(SNR(ii)),'_model',num2str(model_1),'_from_',num2str(model_2),'.mat'];
    recon4 = load(filename4);
    
    % signal 
    est_signal = recon4.dataToMatch_St_unnormal;
    org_signal = recon1.orig_signal;
        
    % estimated values
    S0_est       = double(squeeze(recon1.S0_map));
    T10_est      = double(squeeze(recon1.T10_map));
    Ktrans_est   = double(squeeze(recon2.Ktrans_map));
    vp_est       = double(squeeze(recon2.vp_map));
    ve_est       = double(squeeze(recon2.ve_map));
    BAT_est      = double(squeeze(recon2.BAT_map));
    
    % true values
    S0_true     = recon3.S0_map;
    T10_true    = recon3.T10_map;
    Ktrans_true = recon3.Ktrans_map;
    vp_true     = recon3.vp_map;
    ve_true     = recon3.ve_map;
    BAT_true    = recon3.BAT_map;
    
    % dimensions 
    dim1 = size(S0_true,1);
    dim2 = size(S0_true,2);
    N = dim1*dim2;
    
    RMSE_fit_mean(ii) = RMSE_GoodnessOfFit(est_signal,org_signal,N,0);
    
    % error
    error_Ktrans    = (Ktrans_est - Ktrans_true)./Ktrans_true*100;
    error_vp        = (vp_est - vp_true)./vp_true*100;
    error_ve        = (ve_est - ve_true)./ve_true*100;
    error_T10       = (T10_est - T10_true)./T10_true*100;
    error_S0        = (S0_est - S0_true)./S0_true*100;
    error_BAT       = (BAT_est - BAT_true)./BAT_true*100;
    
    [NRMSE_Ktrans(ii),RMSE_Ktrans(ii)] = NRMSE_calc(Ktrans_true,Ktrans_est,N,0);
    [NRMSE_vp(ii),RMSE_vp(ii)] = NRMSE_calc(vp_true,vp_est,N,0);
    [NRMSE_ve(ii),RMSE_ve(ii)] = NRMSE_calc(ve_true,ve_est,N,0);
    [NRMSE_T10(ii),RMSE_T10(ii)] = NRMSE_calc(T10_true,T10_est,N,0);
    [NRMSE_BAT(ii),RMSE_BAT(ii)] = NRMSE_calc(BAT_true,BAT_est,N,0);
    [NRMSE_S0(ii),RMSE_S0(ii)] = NRMSE_calc(S0_true,S0_est,N,0);
    
    % filter low signal 
    
    if filter_low_signal == 'y'
        [indx_mm,indx_nn,no_zeros_ii,org_signal_filter_low_signal,est_signal_filter_low_signal] = find_low_signal(org_signal,est_signal,cutoff);
        indx_mm_record{ii} = indx_mm;
        indx_nn_record{ii} = indx_nn;
        no_zeros(ii) = no_zeros_ii;

        if indx_mm ~= 0 
            % true
            Ktrans_true_filter_low_signal   = eliminate_low_signal(indx_mm,indx_nn,Ktrans_true);
            vp_true_filter_low_signal       = eliminate_low_signal(indx_mm,indx_nn,vp_true);
            ve_true_filter_low_signal       = eliminate_low_signal(indx_mm,indx_nn,ve_true);
            T10_true_filter_low_signal      = eliminate_low_signal(indx_mm,indx_nn,T10_true);
            S0_true_filter_low_signal       = eliminate_low_signal(indx_mm,indx_nn,S0_true);
            BAT_true_filter_low_signal      = eliminate_low_signal(indx_mm,indx_nn,BAT_true);
            
            % estimate
            Ktrans_est_filter_low_signal    = eliminate_low_signal(indx_mm,indx_nn,Ktrans_est);
            vp_est_filter_low_signal        = eliminate_low_signal(indx_mm,indx_nn,vp_est);
            ve_est_filter_low_signal        = eliminate_low_signal(indx_mm,indx_nn,ve_est);
            T10_est_filter_low_signal       = eliminate_low_signal(indx_mm,indx_nn,T10_est);
            S0_est_filter_low_signal        = eliminate_low_signal(indx_mm,indx_nn,S0_est);
            BAT_est_filter_low_signal       = eliminate_low_signal(indx_mm,indx_nn,BAT_est);
            
            % error
            error_Ktrans_filter_low_signal  = (Ktrans_est_filter_low_signal - Ktrans_true_filter_low_signal)./Ktrans_true_filter_low_signal*100;
            error_vp_filter_low_signal      = (vp_est_filter_low_signal - vp_true_filter_low_signal)./vp_true_filter_low_signal*100;
            error_ve_filter_low_signal      = (ve_est_filter_low_signal - ve_true_filter_low_signal)./ve_true_filter_low_signal*100;
            error_T10_filter_low_signal     = (T10_est_filter_low_signal - T10_true_filter_low_signal)./T10_true_filter_low_signal*100;
            error_S0_filter_low_signal      = (S0_est_filter_low_signal - S0_true_filter_low_signal)./S0_true_filter_low_signal*100;
            error_BAT_filter_low_signal     = (BAT_est_filter_low_signal - BAT_true_filter_low_signal)./BAT_true_filter_low_signal*100;

            [NRMSE_Ktrans_filter_low_signal(ii),RMSE_Ktrans_filter_low_signal(ii)]  = NRMSE_calc(Ktrans_true_filter_low_signal,Ktrans_est_filter_low_signal,N,no_zeros(ii));
            [NRMSE_vp_filter_low_signal(ii),RMSE_vp_filter_low_signal(ii)]          = NRMSE_calc(vp_true_filter_low_signal,vp_est_filter_low_signal,N,no_zeros(ii));
            [NRMSE_ve_filter_low_signal(ii),RMSE_ve_filter_low_signal(ii)]          = NRMSE_calc(ve_true_filter_low_signal,ve_est_filter_low_signal,N,no_zeros(ii));
            [NRMSE_T10_filter_low_signal(ii),RMSE_T10_filter_low_signal(ii)]        = NRMSE_calc(T10_true_filter_low_signal,T10_est_filter_low_signal,N,no_zeros(ii));
            [NRMSE_BAT_filter_low_signal(ii),RMSE_BAT_filter_low_signal(ii)]        = NRMSE_calc(BAT_true_filter_low_signal,BAT_est_filter_low_signal,N,no_zeros(ii));
            [NRMSE_S0_filter_low_signal(ii),RMSE_S0_filter_low_signal(ii)]          = NRMSE_calc(S0_true_filter_low_signal,S0_est_filter_low_signal,N,no_zeros(ii));
        end
        RMSE_fit_mean_filter_low_signal(ii) = RMSE_GoodnessOfFit(est_signal_filter_low_signal,org_signal_filter_low_signal,N,no_zeros(ii));
    end
        
    if filter_large_error == 'y'
    
        error_param             = (eval(['error_',num2str(param_cutoff)]));
        max(error_param(:))
        
        [indx,indx_cutoff]      = ErrorFilterIndex(error_param,param_error_cutoff);
        
        [Ktrans_true_filtered,Ktrans_true_extra_dict]     = ErrorFilter(Ktrans_true,indx,indx_cutoff);
        [vp_true_filtered,vp_true_extra_dict]             = ErrorFilter(vp_true,indx,indx_cutoff);
        [ve_true_filtered,ve_true_extra_dict]             = ErrorFilter(ve_true,indx,indx_cutoff);
        [T10_true_filtered,T10_true_extra_dict]           = ErrorFilter(T10_true,indx,indx_cutoff);
        [BAT_true_filtered,BAT_true_extra_dict]           = ErrorFilter(BAT_true,indx,indx_cutoff);
        [S0_true_filtered,S0_true_extra_dict]             = ErrorFilter(S0_true,indx,indx_cutoff);
        
        [Ktrans_est_filtered,Ktrans_est_extra_dict]       = ErrorFilter(Ktrans_est,indx,indx_cutoff);
        [vp_est_filtered,vp_est_extra_dict]               = ErrorFilter(vp_est,indx,indx_cutoff);
        [ve_est_filtered,ve_est_extra_dict]               = ErrorFilter(ve_est,indx,indx_cutoff);
        [T10_est_filtered,T10_est_extra_dict]             = ErrorFilter(T10_est,indx,indx_cutoff);
        [BAT_est_filtered,BAT_est_extra_dict]             = ErrorFilter(BAT_est,indx,indx_cutoff);
        [S0_est_filtered,S0_est_extra_dict]               = ErrorFilter(S0_est,indx,indx_cutoff);
        
        
%         [Ktrans_true_filtered,Ktrans_est_filtered,Ktrans_error_filtered,param_indx,param_indx_cutoff] = ErrorFilter(Ktrans_true,Ktrans_est,error_param, SNR(ii),param_error_cutoff);
%         [vp_true_filtered,vp_est_filtered,vp_error_filtered,param_indx,param_indx_cutoff] = ErrorFilter(vp_true,vp_est,error_param, SNR(ii),param_error_cutoff);
%         [ve_true_filtered,ve_est_filtered,ve_error_filtered,param_indx,param_indx_cutoff] = ErrorFilter(ve_true,ve_est,error_param, SNR(ii),param_error_cutoff);
%         [T10_true_filtered,T10_est_filtered,T10_error_filtered,param_indx,param_indx_cutoff] = ErrorFilter(T10_true,T10_est,error_param, SNR(ii),param_error_cutoff);
%         [BAT_true_filtered,BAT_est_filtered,BAT_error_filtered,param_indx,param_indx_cutoff] = ErrorFilter(BAT_true,BAT_est,error_param, SNR(ii),param_error_cutoff);
%         [S0_true_filtered,S0_est_filtered,S0_error_filtered,param_indx,param_indx_cutoff] = ErrorFilter(S0_true,S0_est,error_param, SNR(ii),param_error_cutoff);  
        
        N_filtered(ii) = indx_cutoff;
        
        [NRMSE_Ktrans_filtered(ii),RMSE_Ktrans_filtered(ii)]    = NRMSE_calc(Ktrans_true_filtered,Ktrans_est_filtered,N_filtered(ii),0);
        [NRMSE_vp_filtered(ii),RMSE_vp_filtered(ii)]            = NRMSE_calc(vp_true_filtered,vp_est_filtered,N_filtered(ii),0);
        [NRMSE_ve_filtered(ii),RMSE_ve_filtered(ii)]            = NRMSE_calc(ve_true_filtered,ve_est_filtered,N_filtered(ii),0);
        [NRMSE_T10_filtered(ii),RMSE_T10_filtered(ii)]          = NRMSE_calc(T10_true_filtered,T10_est_filtered,N_filtered(ii),0);
        [NRMSE_BAT_filtered(ii),RMSE_BAT_filtered(ii)]          = NRMSE_calc(BAT_true_filtered,BAT_est_filtered,N_filtered(ii),0);
        [NRMSE_S0_filtered(ii),RMSE_S0_filtered(ii)]            = NRMSE_calc(S0_true_filtered,S0_est_filtered,N_filtered(ii),0);
        
        est_signal_filtered = reshape(est_signal,[size(est_signal,1),size(est_signal,2)*size(est_signal,3)]);
        est_signal_filtered = est_signal_filtered(:,indx);
        est_signal_filtered = est_signal_filtered(:,1:indx_cutoff);
        
        org_signal_filtered = reshape(org_signal,[size(org_signal,1),size(org_signal,2)*size(org_signal,3)]);
        org_signal_filtered = org_signal_filtered(:,indx);
        org_signal_filtered = org_signal_filtered(:,1:indx_cutoff);
        
        RMSE_fit_filtred_signal_mean(ii) = RMSE_GoodnessOfFit(est_signal_filtered,org_signal_filtered,indx_cutoff,0);
    end
       
    if plot_all_figures == 'y'
        fig = figure;
        subplot(3,6,1)
        imagesc(Ktrans_true)
        colorbar 
        colormap jet
        axis off
        axis square
        caxis(ktransCaxis)
        title('True K^{trans} map')

        subplot(3,6,7)
        imagesc(Ktrans_est)
        colorbar 
        colormap jet
        axis off
        axis square
        caxis(ktransCaxis)
        title('Est Ktrans map')

        subplot(3,6,13)
        imagesc(squeeze(error_Ktrans))
        colorbar 
        colormap jet
        axis off
        axis square
        caxis(Error_max)
        title('Ktrans Error %')

        subplot(3,6,2)
        imagesc(vp_true)
        colorbar 
        colormap jet
        axis off
        axis square
        caxis(vpCaxis)
        title('True vp map')

        subplot(3,6,8)
        imagesc(vp_est)
        colorbar 
        colormap jet
        axis off
        axis square
        caxis(vpCaxis)
        title('Est vp map')

        subplot(3,6,14)
        imagesc(squeeze(error_vp))
        colorbar 
        colormap jet
        axis off
        axis square
        caxis(Error_max)
        title('vp Error %')

        subplot(3,6,3)
        imagesc(ve_true)
        colorbar 
        colormap jet
        axis off
        axis square
        caxis(veCaxis)
        title('True ve map')

        subplot(3,6,9)
        imagesc(ve_est)
        colorbar 
        colormap jet
        axis off
        axis square
        caxis(veCaxis)
        title('Est ve map')

        subplot(3,6,15)
        imagesc(squeeze(error_ve))
        colorbar 
        colormap jet
        axis off
        axis square
        caxis(Error_max)
        title('ve Error %')

        subplot(3,6,4)
        imagesc(T10_true)
        colorbar 
        colormap jet
        axis off
        axis square
        caxis(T10Caxis)
        title('True T1 map')

        subplot(3,6,10)
        imagesc(T10_est)
        colorbar 
        colormap jet
        axis off
        axis square
        caxis(T10Caxis)
        title('Est T1 map')

        subplot(3,6,16)
        imagesc(squeeze(error_T10))
        colorbar 
        colormap jet
        axis off
        axis square
        caxis(Error_max)
        title('T10 Error %')

        subplot(3,6,5)
        imagesc(S0_true)
        colorbar 
        colormap jet
        axis off
        axis square
        caxis(S0Caxis)
        title('True S0 map')

        subplot(3,6,11)
        imagesc(S0_est)
        colorbar 
        colormap jet
        axis off
        axis square
        caxis(S0Caxis)
        title('Est S0 map')

        subplot(3,6,17)
        imagesc(squeeze(error_S0))
        colorbar 
        colormap jet
        axis off
        axis square
        caxis(Error_max)
        title('S0 Error %')

        subplot(3,6,6)
        imagesc(BAT_true)
        colorbar 
        colormap jet
        axis off
        axis square
        caxis(BATCaxis)
        title('True BAT map')

        subplot(3,6,12)
        imagesc(BAT_est)
        colorbar 
        colormap jet
        axis off
        axis square
        caxis(BATCaxis)
        title('Est BAT map')

        subplot(3,6,18)
        imagesc(squeeze(error_BAT))
        colorbar 
        colormap jet
        axis off
        axis square
        caxis(Error_max)
        title('BAT Error %')

        set(fig,'Units', 'normalized', 'Position', [0.2,0.37,.55,.45])
        sgtitle(['model: ',num2str(model_1),', from model: ',num2str(model_2),', SNR: ',num2str(SNR(ii)),', version',num2str(version)])
    end
%% RMSE and NRMSE calculation for each estimated parameter 
    
    if filter_low_T1 == 'y'
        [T10_map_filtered,est_map_filtered,indx] = filter_low_values(T10_true,T10_est,T1_treshold);
        
        Ktrans_true = Ktrans_true(:);
        vp_true     = vp_true(:);
        ve_true     = ve_true(:);
        T10_true    = T10_true(:);
        BAT_true    = BAT_true(:);
        S0_true     = S0_true(:);
        
        Ktrans_true = Ktrans_true(indx);
        T10_true    = T10_true(indx);
        vp_true     = vp_true(indx);
        ve_true     = ve_true(indx);
        BAT_true    = BAT_true(indx);
        S0_true     = S0_true(indx);
        
        Ktrans_est  = Ktrans_est(indx);
        T10_est     = T10_est(indx);
        vp_est      = vp_est(indx);
        ve_est      = ve_est(indx);
        BAT_est     = BAT_est(indx);
        S0_est      = S0_est(indx);
        N           = length(indx);
    end    
    
%% SSIM

% [ssimval_Ktrans,ssimmap_Ktrans] = ssim(Ktrans_est,Ktrans_true);
% [ssimval_vp,ssimmap_vp] = ssim(vp_est,vp_true);
% [ssimval_ve,ssimmap_ve] = ssim(ve_est,ve_true);
% [ssimval_T10,ssimmap_T10] = ssim(T10_est,T10_true);
% [ssimval_S0,ssimmap_S0] = ssim(S0_est,S0_true);
% [ssimval_BAT,ssimmap_BAT] = ssim(BAT_est,BAT_true);
% 
% %%
% 
% fig = figure,
% subplot(2,3,1)
% imagesc(ssimmap_Ktrans)
% title(['Ktrans SSIM: ',num2str(ssimval_Ktrans)])
% axis off
% colorbar
% set(gca,'FontSize',14)
% 
% subplot(2,3,2)
% imagesc(ssimmap_vp)
% title(['vp SSIM: ',num2str(ssimval_vp)])
% axis off
% colorbar
% set(gca,'FontSize',14)
% 
% subplot(2,3,3)
% imagesc(ssimmap_ve)
% title(['ve SSIM: ',num2str(ssimval_ve)])
% axis off
% colorbar
% set(gca,'FontSize',14)
% 
% subplot(2,3,4)
% imagesc(ssimmap_T10)
% title(['T10 SSIM: ',num2str(ssimval_T10)])
% axis off
% colorbar
% set(gca,'FontSize',14)
% 
% subplot(2,3,5)
% imagesc(ssimmap_S0)
% title(['S0 SSIM: ',num2str(ssimval_S0)])
% axis off
% colorbar
% set(gca,'FontSize',14)
% 
% subplot(2,3,6)
% imagesc(ssimmap_BAT)
% title(['BAT SSIM: ',num2str(ssimval_BAT)])
% axis off
% colorbar
% set(gca,'FontSize',14)
% 
% set(fig,'Units', 'normalized', 'Position', [0.2,0.37,.55,.45])
% sgtitle(['SSIM: AIF Type: ',num2str(AIF_type),', model: ',num2str(model_1),', from model: ',num2str(model_2),', SNR: ',num2str(SNR(ii)),', version',num2str(version)])

%% RMSE and NRMSE calculation for fitted signal vs original signal 
    
    if plot_all_figures == 'y'
    m = 0;
    n = 0;
    while min(m(:)) == 0 || min(n(:)) == 0
        m = round(rand(1,4)*size(est_signal,2));
        n = round(rand(1,4)*size(est_signal,3));
    end
        fig = figure;
        for jj=1:length(m)
            subplot(1,4,jj)
            plot(org_signal(:,m(jj),n(jj)),'k-','linewidth',1.5), hold on
            plot(est_signal(:,m(jj),n(jj)),'r','linewidth',1.5), hold on
            legend('Original signal','Estimated signal')
            title(['SNR: ',num2str(SNR(ii))])
        end
    end
    if filter_large_error == 'y'
        if SNR(ii) == SNR_extra_dict
            
            fig = figure;
            subplot(2,3,1)
            plot(Ktrans_true_extra_dict,Ktrans_est_extra_dict,'.'),title('Ktrans'),
            xlim([0 max(max(Ktrans_true_extra_dict),max(Ktrans_est_extra_dict))])
            ylim([0 max(max(Ktrans_true_extra_dict),max(Ktrans_est_extra_dict))])
            axis square
            xlabel('true')
            ylabel('estimate')
            subplot(2,3,2)
            plot(vp_true_extra_dict,vp_est_extra_dict,'.'),title('vp'),
            xlim([0 max(max(vp_true_extra_dict),max(vp_est_extra_dict))])
            ylim([0 max(max(vp_true_extra_dict),max(vp_est_extra_dict))])
            axis square
            xlabel('true')
            ylabel('estimate')
            subplot(2,3,3)
            plot(ve_true_extra_dict,ve_est_extra_dict,'.'),title('ve'),
            xlim([0 max(max(ve_true_extra_dict),max(ve_est_extra_dict))])
            ylim([0 max(max(ve_true_extra_dict),max(ve_est_extra_dict))])
            axis square
            xlabel('true')
            ylabel('estimate')
            subplot(2,3,4)
            plot(T10_true_extra_dict,T10_est_extra_dict,'.'),title('T10'),
            xlim([0 max(max(T10_true_extra_dict),max(T10_est_extra_dict))])
            ylim([0 max(max(T10_true_extra_dict),max(T10_est_extra_dict))])
            axis square
            xlabel('true')
            ylabel('estimate')
            subplot(2,3,5)
            plot(BAT_true_extra_dict,BAT_est_extra_dict,'.'),title('BAT'),
            xlim([0 max(max(BAT_true_extra_dict),max(BAT_est_extra_dict))])
            ylim([0 max(max(BAT_true_extra_dict),max(BAT_est_extra_dict))])
            axis square
            xlabel('true')
            ylabel('estimate')
            subplot(2,3,6)
            plot(S0_true_extra_dict,S0_est_extra_dict,'.'),title('S0'),
            xlim([0 max(max(S0_true_extra_dict),max(S0_est_extra_dict))])
            ylim([0 max(max(S0_true_extra_dict),max(S0_est_extra_dict))])
            axis square
            xlabel('true')
            ylabel('estimate')

            Ktrans = Ktrans_true_extra_dict;
            vp = vp_true_extra_dict;
            ve = ve_true_extra_dict;
            T10 = T10_true_extra_dict;
            BAT = BAT_true_extra_dict;
            S0 = S0_true_extra_dict;
            B1 = ones(size(S0));
            filename = ['Extra_dict_V_',num2str(version),'_AIF_',num2str(AIF_type),'_model_',num2str(model_1),'_',num2str(model_2),'_LargeErrorCutoff_',num2str(filter_large_error),'_',num2str(param_error_cutoff),'_SNR_',num2str(SNR_extra_dict),'.mat']

            save(filename,'Ktrans','vp','ve','T10','S0','BAT','B1')
        end
    end
end

%% Bar plot for Goodness of fit

bar_width = 0.2;

fig = figure;
subplot(1,3,1)
bar((SNR),RMSE_fit_mean(1:end),bar_width); hold on 
xlabel('SNR')
ylabel('MSE for fitting')
xlim([0 110])
str0 = {'all data'};
legend(str0)
title('nRMSE Goodness of fit')
sgtitle({['V.',num2str(version),' - ',num2str(dictionary_size),' - model: ',num2str(model_1),'/',num2str(model_2),' - AIFType: ',num2str(AIF_type)]})
set(gca,'FontSize',14)

if filter_large_error == 'y'
    subplot(1,3,1)
    bar((SNR+bar_width*5),RMSE_fit_filtred_signal_mean(1:end),bar_width); 
    xlabel('SNR')
    ylabel('MSE for fitting')
    xlim([0 110])
    str0 = cat(2,str0,['filtered data with ',num2str(param_cutoff),' error cutoff: ',num2str(param_error_cutoff),'%']);
    legend(str0)
    set(gca,'FontSize',14)
end

if filter_low_signal == 'y'
    subplot(1,3,1)
    bar((SNR+bar_width*10),RMSE_fit_mean_filter_low_signal(1:end),bar_width);
    xlabel('SNR')
    ylabel('MSE for fitting')
    xlim([0 110])
    str0 = cat(2,str0,['filtered low signal cutoff: ',num2str(cutoff),'%']);
    legend(str0)
    set(gca,'FontSize',14)
end

set(fig,'Units', 'normalized', 'Position', [0.1,0.3,.8,.45])

% plot NRMSE vs. SNR for each parameter

% fig = figure;
subplot(1,3,2:3)
plot(SNR,NRMSE_Ktrans,'o-','color',[0 0.4470 0.7410],'markerfacecolor',[0 0.4470 0.7410],'linewidth',1.5), hold on 
plot(SNR,NRMSE_vp,'*-','color',[0.8500 0.3250 0.0980],'markerfacecolor',[0.8500 0.3250 0.0980],'linewidth',1.5), hold on 
plot(SNR,NRMSE_ve,'<-','color',[0.9290 0.6940 0.1250],'markerfacecolor',[0.9290 0.6940 0.1250]	,'linewidth',1.5), hold on 
plot(SNR,NRMSE_T10,'x-','color',[0.4940 0.1840 0.5560],'markerfacecolor',[0.4940 0.1840 0.5560],'linewidth',1.5), hold on 
plot(SNR,NRMSE_S0,'d-','color',[0.4660 0.6740 0.1880],'markerfacecolor',[0.4660 0.6740 0.1880],'linewidth',1.5), hold on 
plot(SNR,NRMSE_BAT,'^-','color',[0.3010 0.7450 0.9330],'markerfacecolor',[0.3010 0.7450 0.9330],'linewidth',1.5), hold on 
str1 = {'Ktrans','vp','ve','T10','S0','BAT'};
legend(str1,'linewidth',1.5)
xlabel('SNR')
ylabel('nRMSE (%)')
title('nRMSE')
% sgtitle({['V.',num2str(version),' - ',num2str(dictionary_size)];['AIF:',num2str(AIF_type),' - model:',num2str(model_1),'/',num2str(model_2)]})
xlim([0 200])
set(gca,'FontSize',14)

if filter_low_signal == 'y'
    subplot(1,3,2:3)
    plot(SNR,NRMSE_Ktrans_filter_low_signal,'o:','color',[0 0.4470 0.7410],'markerfacecolor',[0 0.4470 0.7410],'linewidth',1.5), hold on 
    plot(SNR,NRMSE_vp_filter_low_signal,'*:','color',[0.8500 0.3250 0.0980],'markerfacecolor',[0.8500 0.3250 0.0980],'linewidth',1.5), hold on 
    plot(SNR,NRMSE_ve_filter_low_signal,'<:','color',[0.9290 0.6940 0.1250],'markerfacecolor',[0.9290 0.6940 0.1250]	,'linewidth',1.5), hold on 
    plot(SNR,NRMSE_T10_filter_low_signal,'x:','color',[0.4940 0.1840 0.5560],'markerfacecolor',[0.4940 0.1840 0.5560],'linewidth',1.5), hold on 
    plot(SNR,NRMSE_S0_filter_low_signal,'d:','color',[0.4660 0.6740 0.1880],'markerfacecolor',[0.4660 0.6740 0.1880],'linewidth',1.5), hold on 
    plot(SNR,NRMSE_BAT_filter_low_signal,'^:','color',[0.3010 0.7450 0.9330],'markerfacecolor',[0.3010 0.7450 0.9330],'linewidth',1.5), hold on 
    str{1} = ['Ktrans for filtered signal - cutoff: ',num2str(cutoff)];
    str{2} = ['vp for filtered signal - cutoff: ',num2str(cutoff)];
    str{3} = ['ve for filtered signal - cutoff: ',num2str(cutoff)];
    str{4} = ['T10 for filtered signal - cutoff: ',num2str(cutoff)];
    str{5} = ['BAT for filtered signal - cutoff: ',num2str(cutoff)];
    str{6} = ['S0 for filtered signal - cutoff: ',num2str(cutoff)];
    str1 = cat(2,str1,str);
    legend(str1)
    set(gca,'FontSize',14)
end

if filter_large_error == 'y'
    subplot(1,3,2:3)
    plot(SNR,NRMSE_Ktrans_filtered,'o--','color',[0 0.4470 0.7410],'markerfacecolor',[0 0.4470 0.7410],'linewidth',1.5), hold on 
    plot(SNR,NRMSE_vp_filtered,'*--','color',[0.8500 0.3250 0.0980],'markerfacecolor',[0.8500 0.3250 0.0980],'linewidth',1.5), hold on 
    plot(SNR,NRMSE_ve_filtered,'<--','color',[0.9290 0.6940 0.1250],'markerfacecolor',[0.9290 0.6940 0.1250]	,'linewidth',1.5), hold on 
    plot(SNR,NRMSE_T10_filtered,'x--','color',[0.4940 0.1840 0.5560],'markerfacecolor',[0.4940 0.1840 0.5560],'linewidth',1.5), hold on 
    plot(SNR,NRMSE_S0_filtered,'d--','color',[0.4660 0.6740 0.1880],'markerfacecolor',[0.4660 0.6740 0.1880],'linewidth',1.5), hold on 
    plot(SNR,NRMSE_BAT_filtered,'^--','color',[0.3010 0.7450 0.9330],'markerfacecolor',[0.3010 0.7450 0.9330],'linewidth',1.5), hold on 
    str2{1} = ['Ktrans for ',num2str(param_cutoff),' error cutoff: ',num2str(param_error_cutoff),'%'];
    str2{2} = ['vp with ',num2str(param_cutoff),' error cutoff: ',num2str(param_error_cutoff),'%'];
    str2{3} = ['ve with ',num2str(param_cutoff),' error cutoff: ',num2str(param_error_cutoff),'%'];
    str2{4} = ['T10 with ',num2str(param_cutoff),' error cutoff: ',num2str(param_error_cutoff),'%'];
    str2{5} = ['BAT with ',num2str(param_cutoff),' error cutoff: ',num2str(param_error_cutoff),'%'];
    str2{6} = ['S0 with ',num2str(param_cutoff),' error cutoff: ',num2str(param_error_cutoff),'%'];
    str1 = cat(2,str1,str2);
    legend(str1)
    set(gca,'FontSize',14)
end





    
    
%%

% filename5 = ['NRMSE_analysis_Version_',num2str(version),'_AIFType_',num2str(AIF_type),'_model_',num2str(model_1),'_from_',num2str(model_2),'_LowSignalCutoff_',num2str(filter_low_signal),'_',num2str(cutoff),'_LargeErrorCutoff_',num2str(filter_large_error),'.mat']
% if filter_large_error == 'n'
%     save(filename5, 'NRMSE_Ktrans','RMSE_Ktrans','NRMSE_vp','RMSE_vp','NRMSE_ve','RMSE_ve','NRMSE_T10','RMSE_T10','NRMSE_S0','RMSE_S0','NRMSE_BAT','RMSE_BAT','SNR','RMSE_fit_mean')
% elseif filter_large_error == 'y'
%     save(filename5, 'NRMSE_Ktrans','RMSE_Ktrans','NRMSE_vp','RMSE_vp','NRMSE_ve','RMSE_ve','NRMSE_T10','RMSE_T10','NRMSE_S0','RMSE_S0','NRMSE_BAT','RMSE_BAT','SNR','RMSE_fit_mean',...
%                     'NRMSE_Ktrans_filtered','RMSE_Ktrans_filtered','NRMSE_vp_filtered','RMSE_vp_filtered','NRMSE_ve_filtered','RMSE_ve_filtered',...
%                     'NRMSE_T10_filtered','RMSE_T10_filtered','NRMSE_BAT_filtered','RMSE_BAT_filtered','NRMSE_S0_filtered','RMSE_S0_filtered','RMSE_fit_filtred_signal_mean')
% end    
% 
% %%
% 
% plot(SNR,NRMSE_Ktrans,'o--','color',[0 0.4470 0.7410],'markerfacecolor',[0 0.4470 0.7410],'linewidth',1.5), hold on 
% plot(SNR,NRMSE_vp,'*--','color',[0.8500 0.3250 0.0980],'markerfacecolor',[0.8500 0.3250 0.0980],'linewidth',1.5), hold on 
% plot(SNR,NRMSE_ve,'<--','color',[0.9290 0.6940 0.1250],'markerfacecolor',[0.9290 0.6940 0.1250]	,'linewidth',1.5), hold on 
% plot(SNR,NRMSE_T10,'x--','color',[0.4940 0.1840 0.5560],'markerfacecolor',[0.4940 0.1840 0.5560],'linewidth',1.5), hold on 
% plot(SNR,NRMSE_S0,'d--','color',[0.4660 0.6740 0.1880],'markerfacecolor',[0.4660 0.6740 0.1880],'linewidth',1.5), hold on 
% plot(SNR,NRMSE_BAT,'^--','color',[0.3010 0.7450 0.9330],'markerfacecolor',[0.3010 0.7450 0.9330],'linewidth',1.5), hold on 
% %%
% legend('Ktrans','vp','ve','T10','S0','BAT','linewidth',1.5)
% xlabel('SNR')
% ylabel('nRMSE (%)')
% title(['nRMSE Analysis Numeric Sim V.',num2str(version),' - AIFType:',num2str(AIF_type),' - model:',num2str(model_1),'/',num2str(model_2),' - Cutoff: ',num2str(filter_low_signal),' ',num2str(cutoff)])
% 
% xlim([0 200])
% set(gca,'FontSize',14)
% 
% %%
% % fig = figure
% % % boxplot(RMSE_fit_all,'Notch','on','Labels',{'1000','100'})
% % boxplot(RMSE_fit_all,'Notch','on','Labels',{'1000','200','100','90','80','70','60','50','40','30','20','10','5'})
% % xlabel('SNR')
% % ylim([0 1])
% % title(['nRMSE Analysis Goodness of fit Numeric Sim V.',num2str(version),' - AIFType:',num2str(AIF_type),' - model:',num2str(model_1),'/',num2str(model_2)])
% % set(gca,'FontSize',14)
% 
% %%
% 
% % T10_true(m(2),n(2))
% % T10_est(m(2),n(2))
% % 
% % Ktrans_true(m(2),n(2))
% % Ktrans_est(m(2),n(2))
% % 
% % vp_true(m(2),n(2))
% % vp_est(m(2),n(2))
% % 
% % ve_true(m(2),n(2))
% % ve_est(m(2),n(2))
% % 
% % S0_true(m(2),n(2))
% % S0_est(m(2),n(2))
% % 
% % BAT_true(m(2),n(2))
% % BAT_est(m(2),n(2))
% % 
% % a1_true(m(2),n(2))
% % a1_est(m(2),n(2))
% 







