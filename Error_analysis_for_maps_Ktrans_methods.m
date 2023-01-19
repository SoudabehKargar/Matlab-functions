clc
clear variables
close all
 
%%

all_SNR = 'n';
if all_SNR == 'n'
    SNR = [1000]'
else
    SNR = [1000,200,100,90,80,70,60,50,40,30,20,10]'
end

dictionary_size = '200k'
AIF = 'P'              %  P: Parker M: Modified Parker S: StepLike

if dictionary_size == '400k'    
    if AIF == 'S'      % 282/280  289/290
        version = 51
        model_1 = 290
        model_2 = 289

    elseif AIF == 'P'  % 254/256 291/292
        version = 50
        model_1 = 292
        model_2 = 291

    elseif AIF == 'M'  % 275/281 287/288
        version = 52
        model_1 = 288
        model_2 = 287
    end
elseif dictionary_size == '200k'
    if AIF == 'S'      % 282/280  289/290
        version = 51
        model_1 = 296
        model_2 = 295

    elseif AIF == 'P'  % 254/256 291/292
        version = 53
        model_1 = 298
        model_2 = 297

    elseif AIF == 'M'  %  275/281 287/288
        version = 54
        model_1 = 294
        model_2 = 293
    end
end 

%% True values
filename3 = ['/Volumes/MRIClinical/kargar/DL/DCE_DRONE_code/test_data_vals/test_data_vals_v',num2str(version),'.mat'];
recon3 = load(filename3);
Ktrans_true  = recon3.Ktrans_map;
vp_true      = recon3.vp_map;
ve_true      = recon3.ve_map;
T10_true     = recon3.T10_map;
S0_true      = recon3.S0_map;
BAT_true     = (recon3.BAT_map);

%% Estimated values 

for ii=1:length(SNR) 
    filename1 = ['/Volumes/MRIClinical/kargar/DL/DCE_DRONE_code/test_data_recon/test_data_recon_num_estimate_test_vals_v',num2str(version),'_SNR',num2str(SNR(ii)),'_model_',num2str(model_2),'_from_.mat'];
    recon1 = load(filename1);
    filename2 = ['/Volumes/MRIClinical/kargar/DL/DCE_DRONE_code/test_data_recon/test_data_recon_num_estimate_test_vals_v',num2str(version),'_SNR',num2str(SNR(ii)),'_model_',num2str(model_1),'_from_',num2str(model_2),'.mat'];
    recon2 = load(filename2);
    
    S0_est(ii,:,:)          = double(squeeze(recon1.S0_map));
    T10_est(ii,:,:)         = double(squeeze(recon1.T10_map));
    Ktrans_est(ii,:,:)      = double(squeeze(recon2.Ktrans_map));
    vp_est(ii,:,:)          = double(squeeze(recon2.vp_map));
    ve_est(ii,:,:)          = double(squeeze(recon2.ve_map));
    BAT_est(ii,:,:)         = double(squeeze(recon2.BAT_map));
    
    error_Ktrans(ii,:,:)    = (squeeze(Ktrans_est(ii,:,:)) - Ktrans_true)./Ktrans_true*100;
    error_vp(ii,:,:)        = (squeeze(vp_est(ii,:,:)) - vp_true)./vp_true*100;
    error_ve(ii,:,:)        = (squeeze(ve_est(ii,:,:)) - ve_true)./ve_true*100;
    error_T10(ii,:,:)       = (squeeze(T10_est(ii,:,:)) - T10_true)./T10_true*100;
    error_S0(ii,:,:)        = (squeeze(S0_est(ii,:,:)) - S0_true)./S0_true*100;
    error_BAT(ii,:,:)       = (squeeze(BAT_est(ii,:,:)) - BAT_true)./BAT_true*100;
end

%% Plot limits
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
BATErrorCaxis = [0 ceil(max(error_BAT(:)))];
S0ErrorCaxis = [0 ceil(max(error_S0(:)))];

Ktrans_error_bin_width = 10;
vp_error_bin_width = 10;
ve_error_bin_width = 10;
T10_error_bin_width = 5;
BAT_error_bin_width = ceil(max(error_BAT(:)))/6;
S0_error_bin_width = ceil(max(error_S0(:)))/6;

Ktrans_error_bin_max = 60;
vp_error_bin_max = 60;
ve_error_bin_max = 60;
T10_error_bin_max = 30;
BAT_error_bin_max = .7*max(error_BAT(:));
S0_error_bin_max = .7*max(error_S0(:));

%% 
Ktrans_bin_count = ErrorSort(Ktrans_true,  Ktrans_est, error_Ktrans,   SNR,    KtransCaxis,    KtransErrorCaxis,   'Ktrans',   [Ktrans_error_bin_width:    Ktrans_error_bin_width: Ktrans_error_bin_max])
vp_bin_count     = ErrorSort(vp_true,      vp_est,     error_vp,       SNR,    vpCaxis,        vpErrorCaxis,       'vp',       [vp_error_bin_width:        vp_error_bin_width:     vp_error_bin_max])
ve_bin_count     = ErrorSort(ve_true,      ve_est,     error_ve,       SNR,    veCaxis,        veErrorCaxis,       've',       [ve_error_bin_width:        ve_error_bin_width:     ve_error_bin_max])
T10_bin_count    = ErrorSort(T10_true,     T10_est,    error_T10,      SNR,    T10Caxis,       T10ErrorCaxis,      'T10',      [T10_error_bin_width:       T10_error_bin_width:    T10_error_bin_max])
BAT_bin_count    = ErrorSort(BAT_true,     BAT_est,    error_BAT,      SNR,    BATCaxis,       BATErrorCaxis,      'BAT',      [BAT_error_bin_width:       BAT_error_bin_width:    BAT_error_bin_max])
S0_bin_count     = ErrorSort(S0_true,      S0_est,     error_S0,       SNR,    S0Caxis,        S0ErrorCaxis,       'S0',       [S0_error_bin_width:        S0_error_bin_width:     S0_error_bin_max])

%%










