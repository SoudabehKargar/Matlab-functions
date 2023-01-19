function RMSE_fit_mean = RMSE_GoodnessOfFit(est_signal,org_signal,N,no_zeros)

    RMSE_fit = (est_signal - org_signal);
    RMSE_fit = reshape(RMSE_fit,[size(RMSE_fit,1),size(RMSE_fit,2)*size(RMSE_fit,3)]);
    RMSE_fit = RMSE_fit.^2;
    RMSE_fit = sum(RMSE_fit,1);
    RMSE_fit_mean = (sum(RMSE_fit))/(N-no_zeros);
    RMSE_fit_all = RMSE_fit';
    
end
