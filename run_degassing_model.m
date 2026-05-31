function output = run_degassing_model(p)
   

   par_solver_new;
   [CC_growth_input,Thermal_input,H2O_input,Xe_input,Simulation_input,FIG,FD_MOR,FD_p] = input_degassing(p);


    output_degassing = ...
        Degassing_model(Simulation_input,...
                     Thermal_input,CC_growth_input,H2O_input,Xe_input,FD_MOR,FD_p);

   
    %% RMSE

    % observation of degassing rate
    P_mor_pd = p.P_mor_pd;
    Pmor_center = (P_mor_pd(1)+P_mor_pd(2))/2;
    Pmor_range = abs((P_mor_pd(1)-P_mor_pd(2))/2);
    Er_Pmor_pd = (output_degassing.Pmor(end)-Pmor_center)/Pmor_range;



    % observation of Xe: present day concentration & isotope ratio
    if output_degassing.SwXe

        % observation of Xe 1: present day mantle Xe concentration
        center = sum(p.Xe_obs.Xe)/2;
        range = abs(p.Xe_obs.Xe(2)-p.Xe_obs.Xe(1))/2;
        result = output_degassing.Xe(end);
        % RMSE_xec = sqrt( ( (result-center).^2) / (range.^2) );
        Er_xec =   (result-center) ./ range;

        % observation of Xe 2: present day mantle 128Xe/130Xe
        center = sum(p.Xe_obs.Xe128v130)/2;
        range = abs(p.Xe_obs.Xe128v130(2)-p.Xe_obs.Xe128v130(1))/2;
        result = output_degassing.Xe_Atm(end)./ output_degassing.Xe(end);
        % RMSE_xe128v130 = sqrt( ( (result-center).^2) / (range.^2) );
        Er_xe128v130 =   (result-center) ./ range;

        % observation of Xe 3: present day mantle 128Xe/132Xe
        center = sum(p.Xe_obs.Xe128v132)/2;
        range = abs(p.Xe_obs.Xe128v132(2)-p.Xe_obs.Xe128v132(1))/2;
        result = output_degassing.Xe_Atm(end)./ output_degassing.Xe132(end);
        % RMSE_xe128v132 = sqrt( ( (result-center).^2) / (range.^2) );
        Er_xe128v132 =   (result-center) ./ range;
        
        % observation of Xe 4: present day mantle 130Xe/132Xe
        center = sum(p.Xe_obs.Xe130v132)/2;
        range = abs(p.Xe_obs.Xe130v132(2)-p.Xe_obs.Xe130v132(1))/2;
        result = output_degassing.Xe(end)./ output_degassing.Xe132(end);
        % RMSE_xe130v132 = sqrt( ( (result-center).^2) / (range.^2) );
        Er_xe130v132 =   (result-center) ./ range;

        % observation of Xe 5: present day mantle 131Xe/132Xe
        center = sum(p.Xe_obs.Xe131v132)/2;
        range = abs(p.Xe_obs.Xe131v132(2)-p.Xe_obs.Xe131v132(1))/2;
        result = output_degassing.Xe131(end)./ output_degassing.Xe132(end);
        %RMSE_xe131v132 = sqrt( ( (result-center).^2) / (range.^2) );
        Er_xe131v132 =   (result-center) ./ range;

        % observation of Xe 6: present day mantle 134Xe/132Xe
        center = sum(p.Xe_obs.Xe134v132)/2;
        range = abs(p.Xe_obs.Xe134v132(2)-p.Xe_obs.Xe134v132(1))/2;
        result = output_degassing.Xe134(end)./ output_degassing.Xe132(end);
        %RMSE_xe134v132 = sqrt( ( (result-center).^2) / (range.^2) );
        Er_xe134v132 =   (result-center) ./ range;

        % observation of Xe 7: present day mantle 136Xe/132Xe
        center = sum(p.Xe_obs.Xe136v132)/2;
        range = abs(p.Xe_obs.Xe136v132(2)-p.Xe_obs.Xe136v132(1))/2;
        result = output_degassing.Xe136(end)./ output_degassing.Xe132(end);
        %RMSE_xe136v132 = sqrt( ( (result-center).^2) / (range.^2) );
        Er_xe136v132 =   (result-center) ./ range;

    else
        % observation of Xe 1: present day mantle Xe concentration
        RMSE_xec = 0;
        Er_xec = 0;
        % observation of Xe 2: present day mantle 128Xe/130Xe
        RMSE_xe128v130 = 0;
        Er_xe128v130 = 0;
        % observation of Xe 3: present day mantle 128Xe/130Xe
        RMSE_xe128v132 = 0;
        Er_xe128v132 = 0;
        % observation of Xe 4: present day mantle 130Xe/132Xe
        RMSE_xe130v132 = 0;
        Er_xe130v132 = 0;
        % observation of Xe 5: present day mantle 131Xe/132Xe
        RMSE_xe131v132 = 0;
        Er_xe131v132 = 0;
        % observation of Xe 6: present day mantle 134Xe/132Xe
        RMSE_xe134v132 = 0;
        Er_xe134v132 = 0;
        % observation of Xe 7: present day mantle 136Xe/132Xe
        RMSE_xe136v132 = 0;
        Er_xe136v132 = 0;
    end


%     output = struct('kr', p.kappa_r, 'kg', p.kappa_g, ...
%                     'Rs', p.Rs, 'Rp', p.Rp, 'ts', p.ts, 'RMSE', RMSE);

    output = [p.time_cc_input,p.kappa_gcc,p.kappa_rcc,p.Rs,p.Rp,... % ===== 5 parameters for CC model =====
        p.Ubse,p.Qc_pd_input,p.dQc_input,p.Ti,p.eta_ref,... % ===== 5 parameters for thermal model =====
        p.Xm_init,p.Ms_init,p.Fr_w,... % ===== 3 parameters for water model =====
        p.Xe_init,p.Fr_Xe,... % ===== 2 parameters for Xe model =====
        p.Fd_mor,p.Fd_p,... % ===== 2 parameters for MOR & Plume degassing model =====
        Er_Pmor_pd,...  % ===== 1 rmse for present day mantle processing rate ===== % 0,...% Er_water1,...  % ===== 1 rmse for water model =====
        Er_xec,Er_xe128v130,Er_xe128v132,Er_xe130v132,Er_xe131v132,Er_xe134v132,Er_xe136v132];
        % RMSE_xec,RMSE_xe128v130,RMSE_xe128v132,RMSE_xe130v132,RMSE_xe131v132,RMSE_xe134v132,RMSE_xe136v132]; 
       % ===== 7 rmse for Xe model =====
end