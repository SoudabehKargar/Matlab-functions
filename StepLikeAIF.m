% clc
% close all
% clear variables
% 
% time_series_length = 60
% Delta_t = 5
% amplitude = 5
% BAT = 10
% Cp = StepLike_AIF(time_series_length,Delta_t,amplitude)
% Cp = apply_bolus_arrival_time_delay(Cp,BAT)

function Cp = StepLike_AIF(time_series_length,Delta_t,amplitude)
    
    % converted from the Python Code for DRONE-DCE **kargar
    tt = linspace(1,time_series_length,time_series_length)-1;
    tsec = tt*Delta_t;
    t_min = tsec/60; % minutes
    Hct = 0.42; % hematocrit
    %% step like AIF
    
    Cb = amplitude*ones(size(t_min));
    
    % according to Parker et al to get the actual Cp we need to divide by the 
    % hematocrit but the matlab code doesn't do that for some reason 
    Cp = Cb/(1-Hct);
end
