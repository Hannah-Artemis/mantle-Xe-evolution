function output = run_cc_growth_model(p)
   
    [Mc_model, Mdd_model, Mud_model, Krw_model] = CC_growth_fun2(p.t, p.ts, p.tmax, ...
        p.Mcp, p.kappa_g, p.Rp, p.Rs, p.kappa_r,p.frw);
    
    [F_model, S_model, m_tp,m,Krw_real_model,s] = Formation_surface_age_fun(p.t, Mud_model, Mdd_model, Mc_model, Krw_model);
    
    RMSE_F = sqrt(sum((F_model - p.F_obs).^2) / sum(p.F_obs.^2));
    RMSE_S = sqrt(sum((S_model - p.S_obs).^2) / sum(p.S_obs.^2));

    output = [p.ts,p.kappa_g,p.kappa_r,p.Rs,p.Rp,p.frw, RMSE_F,RMSE_S];

end