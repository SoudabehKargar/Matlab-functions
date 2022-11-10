function Ct = ToftsModel_DRONE(Ktrans_map,ve_map,vp_map,Cp_delayed,Delta_t,time_series_length)

    
    Ktrans_map = Ktrans_map(:);
    Ct = zeros(length(Ktrans_map),time_series_length);
    ve_map = ve_map(:);
    for ii=1:length(Ktrans_map)
        kep = Ktrans_map(ii)./ve_map(ii);
        Ct(ii,:) = vp_map(ii)*Cp_delayed(:,ii) + Ktrans_map(ii) * expConv_kargar(Cp_delayed(:,ii), kep,Delta_t,time_series_length);

    end
end
