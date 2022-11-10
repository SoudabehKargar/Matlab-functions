function ff = expConv_kargar(A, B,Delta_t,time_series_length)
    % converted from the Python Code for DRONE-DCE **kargar

        % A : Cp or AIF
        % B : kep
        % f : AIF convolved with exp(-kep.t)

        % f = np.zeros(len(A))
        ff = zeros(length(A),1);
        tt =  linspace(1,time_series_length,time_series_length)-1;
        time = (Delta_t/60) *tt; % OC 11/22/2021 convert to minutes
        Ht = exp(-B*time);
        AA = convmatrix(A);
        AA = AA*Delta_t/60;
        % AA = np.transpose(AA)
        for ii = 1:length(A)
            ff(ii) = sum(AA(:,ii).*Ht');
        % f = np.matmul(Ht,AA)
        end
end

        
