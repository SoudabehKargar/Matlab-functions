function Cp = apply_bolus_arrival_time_delay(Cp,BAT)
        %"""
        % apply the bolus arrival time delay to the AIF by shifting the
        % AIF (i.e. Cp) curve. The signal is shifted and zero-padded.
        % """
        BAT = BAT(:);
        shift = round(BAT); % The BAT needs to be integers so we just round it to the nearest one
        shifted_Cp = zeros(length(Cp),size(BAT,1)); %pre-allocate N x P array for each of the P tissue combinations
        for jj = 1:(length(shift))
            if shift(jj) ~= 0
                curr_indx =  shift(jj); % convert to scalar 
                shifted_Cp(curr_indx+1:end,jj) = Cp(1:end-curr_indx);
            else
                shifted_Cp(:,jj) = Cp;
            end
        end
        % write back to Cp
        Cp = shifted_Cp;

end