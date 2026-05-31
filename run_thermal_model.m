function output = run_thermal_model(p)
    
    par_solver_new;

    [Tp,Tm,Qs,Qc,H,Ra,eta] = Thermal_model(p.t,p.Ti,p.eta_ref,p.Qc_pd_input,p.dQc_input,...
    p.BSE,p.Ubse,p.kappa_gcc,p.time_cc_input,p.para_Scale,p.Hsf);


   
    % observation 1: present day potential temperature around 1350C
    RMSE_pd = sqrt( ( (Tp(end)-p.Tp_obs1.Tp_pd).^2) / (p.Tp_obs1.Tp_range.^2) );

    % observation 2: potential temperature history
    i1=round(p.Tp_obs2.t_anchorHerz1/p.dt);
    i2=round(p.Tp_obs2.t_anchorHerz2/p.dt);
    i3=round(p.Tp_obs2.t_anchorHerz3/p.dt);
    i4=round(p.Tp_obs2.t_anchorHerz4/p.dt);
    % % 70,90,100,130 from data range
    RMSE_t = sqrt((((Tp(i1)-p.Tp_obs2.Tp_anchorHerz1)/70)^2 ...
        +((Tp(i2)-p.Tp_obs2.Tp_anchorHerz2)/90)^2 ...
        +((Tp(i3)-p.Tp_obs2.Tp_anchorHerz3)/100)^2 ...
        +((Tp(i4)-p.Tp_obs2.Tp_anchorHerz4)/130)^2)/4);

    % observation 3: present day heat flux within [26.7e12,34e12]
    Qs0 = sum(p.Qs_obs)/2;dQs=abs(p.Qs_obs(2)-p.Qs_obs(1))/2;
    RMSE_Q = sqrt( ( (Qs(end)-Qs0).^2 ) / (dQs.^2) );

%     output = struct('kr', p.kappa_r, 'kg', p.kappa_g, ...
%                     'Rs', p.Rs, 'Rp', p.Rp, 'ts', p.ts, 'RMSE', RMSE);

    output = [p.time_cc_input,p.kappa_gcc,p.kappa_rcc,p.Rs,p.Rp,p.Ubse,p.Qc_pd_input,p.dQc_input,p.Ti,p.eta_ref,RMSE_pd,RMSE_t,RMSE_Q];

end