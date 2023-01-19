function [bin_count] = ErrorSort(param_true,param_est,error_param, SNR,paramCaxis,ErrorCaxis,Title,bins)
    
    % ErrorSort     : This function sorts the errors for each parameter and
    %                 plot the true and estimated parameters based on the bins to
    %                 colorcoded clusters. 
    % input:
    % --------------
    % param_true    : 
    % param_est     : 
    % error_param   :
    % SNR           :
    % paramCaxis    :
    % ErrorCaxis    :
    % Title         :
    % bins          :
    
    % output:
    % --------------
    % bin_count     :

    bins0 = [0, bins];

    for ii=1:length(SNR)
        param_error = error_param(ii,:,:);
        param_error = abs(param_error(:));
        [param_error_sorted,indx] = sort(param_error);

        param_est_ii = param_est(ii,:,:);
        param_est_ii = param_est_ii(:);
        param_est_ii_sorted = param_est_ii(indx);

        param_true = param_true(:);
        param_true_sorted = param_true(indx);

        for jj = 1:length(bins)
            b1(jj) = find(param_error_sorted>bins(jj),1);   % finds the 1st index that has values larger than bins(jj)
        end
        b_end = length(param_true_sorted);
        b = [0 b1 b_end];
        bin_count = diff(b);                                % calculates the number of points in each bin
    end
    fig = figure;
    subplot(2,2,1),plot(param_est_ii,'.'),title(['Estimated ',Title]),ylim(paramCaxis)
    subplot(2,2,2),plot(param_true,'.'),title(['True ',Title]),ylim(paramCaxis)
    subplot(2,2,3),plot(param_est_ii_sorted,'.'),title(['Sorted Estimated ',Title]),ylim(paramCaxis)
    subplot(2,2,4),plot(param_true_sorted,'.'),title(['Sorted True ',Title]),ylim(paramCaxis)
    set(fig,'Units', 'normalized', 'Position', [0.3,0.2,.4,.5])
    %% Plot True and Estimated values vs. Error %

    fig = figure;
    for jj = 1:length(b)-1
        
        subplot(1,2,1)
        plot(squeeze(param_error_sorted(b(jj)+1:b(jj+1))),squeeze(param_est_ii_sorted(b(jj)+1:b(jj+1))),'.'), hold on;
        axis square
        xlim(ErrorCaxis)
        ylim(paramCaxis)
        title(['SNR: ',num2str(SNR(ii))])
        xlabel('Error %')
        ylabel([Title,' (Estimated)'])
        set(fig,'Units', 'normalized', 'Position', [0.2,0.2,.6,.5])
        if jj<length(b)-1
            str{jj} = [num2str(bin_count(jj)),' : ',Title, ' error % : ', num2str(bins0(jj)), ' - ',num2str(bins0(jj+1)) ];
        else
            str{jj} = [num2str(bin_count(jj)),' : ',Title, ' error % > ', num2str(bins0(jj)) ];
        end
        if ii == 1 && jj==length(b)-1
            legend(str)
        end
        set(gca,'FontSize',14)
        
        subplot(1,2,2)
        plot(squeeze(param_error_sorted(b(jj)+1:b(jj+1))),squeeze(param_true_sorted(b(jj)+1:b(jj+1))),'.'), hold on;
        axis square
        xlim(ErrorCaxis)
        ylim(paramCaxis)
        title(['SNR: ',num2str(SNR(ii))])
        xlabel('Error %')
        ylabel([Title, ' (True)'])
        set(fig,'Units', 'normalized', 'Position', [0.2,0.2,.6,.5])
        if jj<length(b)-1
            str1{jj} = [num2str(bin_count(jj)),' : ',Title, ' error % : ', num2str(bins0(jj)), ' - ',num2str(bins0(jj+1)) ];
        else
            str1{jj} = [num2str(bin_count(jj)),' : ',Title, ' error % > ', num2str(bins0(jj)) ];
        end
        if ii == 1 && jj==length(b)-1
            legend(str1)
        end
        set(gca,'FontSize',14)
    end

end





