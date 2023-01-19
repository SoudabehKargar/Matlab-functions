function [NRMSE_param,RMSE_param] = NRMSE_calc(param_true,param_est,N,no_zeros)

    % Ktrans
    m_param = sum(param_true(:))/(N-no_zeros);

%     mean_Ktrans_true = m_Ktrans;
    a_param = (param_est - param_true).^2;
    RMSE_param = sqrt(sum(a_param(:))/(N-no_zeros));
    NRMSE_param = RMSE_param/m_param*100;
    
end
