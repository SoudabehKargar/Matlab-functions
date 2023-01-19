function [param_filtered,param_leftover] = ErrorFilter(param,indx,indx_cutoff)
    
    % ErrorFilter   : This function sorts the errors for each parameter and
    %                 plot the true and estimated parameters based on the bins to
    %                 colorcoded clusters. 
    % input:
    % --------------
    % param         : 
    % indx          : 
    % indx_cutoff   :
    
    % output:
    % --------------
    % param_filtered   :
    
    param           = param(:);
    param_sorted    = param(indx);
    param_filtered  = param_sorted(1:indx_cutoff);
    param_leftover  = param_sorted(indx_cutoff+1:end);
end

   



