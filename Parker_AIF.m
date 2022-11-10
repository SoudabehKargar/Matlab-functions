function Cp = Parker_AIF(time_series_length,Delta_t)
    
    % converted from the Python Code for DRONE-DCE **kargar
    tt = linspace(1,time_series_length,time_series_length)-1;
    tsec = tt*Delta_t;
    t_min = tsec/60; % minutes
    A = [0.809, 0.330]; % scaling const, mmol*min
    T = [0.17046, 0.365]; % centers, min
    sigma = [0.0563, 0.132]; % gaussian width, min
    alpha = 1.050; % amplitude constant of exponential, mmol
    beta = 0.1685; % decay constant of exponential, 1/min
    s = 38.078; % sigmoid width, 1/min
    tau = 0.483; % sigmoid center, min
    Hct = 0.42; % hematocrit
    C = sqrt(2*pi());

    % equation (2) in Parket et al                 
    exp_term0 = ( A(1) / (C*sigma(1)) ) * exp( -(t_min-T(1)) .* (t_min-T(1)) / (2*sigma(1)^2) );
    exp_term1 = ( A(2) / (C*sigma(2)) ) * exp( -(t_min-T(2)) .* (t_min-T(2)) / (2*sigma(2)^2) );
    sig_term = alpha * exp(-beta*t_min) ./ (1 + exp( -s * (t_min-tau) ) );

    Cb = exp_term0 + exp_term1 + sig_term;                       

    % according to Parker et al to get the actual Cp we need to divide by the 
    % hematocrit but the matlab code doesn't do that for some reason 
    Cp = Cb/(1-Hct);
end

