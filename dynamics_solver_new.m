%% outputs values using t, T and Xm


function dynamic_quantity=...
    dynamics_solver_new(t,Xm,T,Thermal_input,CC_growth_input,H2O_input,Xe_input,FD_MOR,FD_p)

    par_solver_new; % read constant inputs including Xp= 2.056e3

    % thermal model related parameters
    T_init = Thermal_input.T_init;
    BSE = Thermal_input.BSE;
    h_El = Thermal_input.h_El;
    Qc_pd = Thermal_input.Qc_pd;
    dQc =Thermal_input.dQc;
    Hsf = Thermal_input.Hsf;
    etain = Thermal_input.etain;
    etalawin =  Thermal_input.etalawin;
    para_Scale = Thermal_input.para_Scale;
    beta = Thermal_input.beta;
    pU = Thermal_input.pU;
    U_cm_yr = Thermal_input.U_cm_yr;
    
    
    
    % CC grwoth model related parametes
    time_cc = CC_growth_input.time_cc;
    kcc = CC_growth_input.kcc;
    krcc = CC_growth_input.krcc;
    Rccs = CC_growth_input.Rccs;
    Rccp = CC_growth_input.Rccp;
    
    % water
    Xm_init = H2O_input.Xm_init;
    Ms_init = H2O_input.Ms_init;
    Fr = H2O_input.Frw;
    Fd_sw = H2O_input.Fd_sw;
    Fd_new = H2O_input.Fd_new;
    f_WI = H2O_input.f_WI;
    
    % Xe
    SwXe = Xe_input.SwXe;
    Xe_init = Xe_input.Xe_init;
    Xes_init = Xe_input.Xes_init;
    XeAtm_init = Xe_input.XeAtm_init;
    XesAtm = Xe_input.XesAtm_init;
    Xe136r130atm_init = Xe_input.Xe136r130atm_init; 
    Xe136_init = Xe_input.Xe136_init;
    Xe131r130atm_init = Xe_input.Xe131r130atm_init;
    Xe131_init = Xe_input.Xe131_init;
    Xe132r130atm_init = Xe_input.Xe132r130atm_init;
    Xe132_init = Xe_input.Xe132_init;
    Xe134r130atm_init = Xe_input.Xe134r130atm_init;
    Xe134_init = Xe_input.Xe134_init;
    Pu_init = Xe_input.Pu_init;
    Ur_init = Xe_input.Ur_init;
    R0Pu_Ur = Xe_input.R0Pu_Ur;
    lamUr = Xe_input.lamUr;
    lamPu = Xe_input.lamPu;
    FrXe = Xe_input.FrXe;
    
    
    Fd=FD_MOR;
    % FD_p = FD_p;

 
    cc = (t<time_cc).*0+(t>=time_cc).*(1-exp(-kcc*(t-time_cc)/yr_s/1e9))./(1-exp(-kcc*(t_pd-time_cc)/yr_s/1e9 ) ) ;

    NG_cc = 1/1e9/yr_s*( (t<time_cc)*0+(t>=time_cc)... %net CC growth rate; Urcc and Urt in atoms/g; X(1) in ppm; thus Xmdot = ppm/s
        .*( kcc*exp(-kcc*(t-time_cc)/yr_s/1e9)...
        ./(1-exp(-kcc*(t_pd-time_cc)/yr_s/1e9 ) ) ) ) .* mcc_val;

    R_cc = 1/1e9/yr_s*( (t<time_cc)*0+(t>=time_cc).* ...
            ( Rccs+ (Rccp-Rccs).*(1-exp(-krcc*(t-time_cc)/yr_s/1e9))./(1-exp(-krcc.*(t_pd-time_cc)/yr_s/1e9)) )); %CC recycling

    Hsf_t=1;
    H_El(:,1) = (h_El(1,1)-U238CC*mcc_val*(cc)/rho/V/Hsf_t).*h_El(2,1).*exp((log(2).*(t_pd-t))./h_El(3,1));
    H_El(:,2) = (h_El(1,2)-U235CC*mcc_val*(cc)/rho/V/Hsf_t).*h_El(2,2).*exp((log(2).*(t_pd-t))./h_El(3,2));
    H_El(:,3) = (h_El(1,3)-ThCC*mcc_val*(cc)/rho/V/Hsf_t).*h_El(2,3).*exp((log(2).*(t_pd-t))./h_El(3,3));
    H_El(:,4) = (h_El(1,4)-KCC*mcc_val*(cc)/rho/V/Hsf_t).*h_El(2,4).*exp((log(2).*(t_pd-t))./h_El(3,4));
   
    Umt=(h_El(1,1)-U238CC*mcc_val*(cc)/rho/V)... % BSE - CC extraction 
        .*exp(lamUr.*(t_pd-t)/yr_s)... % recover decay
        /238*NA; % from weight ratio (U238cc) to atoms/g (Urcc)

    % core heat flux scale to H
    Qc = (Hsf-1)*rho*V .* sum(H_El,2);

    eta_ref = etain;
    eta = eta_ref.*exp(E./(Rg.*T)); % here T is Tavg
    
    if beta == 1/2
        Hsf_rho_V = Hsf*rho*V;
        H  = Hsf_rho_V .* sum(H_El,2); %W %2 in sum() sums along each row not column
        Ra = (alpha_ra.*rhoa.*rhoa.*g.*H.*d^5)./(Kc.*eta.*k);
    else
        Ra = (alpha_ra.*rhoa.*g.*(T-Ts).*d^3)./(k.*eta);        
    end
    
    %U = (k/d).*((L/(pi*d))^(1/3)).*Ra.^(beta); %; 5 .* 0.01 ./ (3600*24*365.25) .* ones(size(T))
    switch para_Scale
        case 1 % Crowley, 2011: using weakly nonlinear theory from Geodynamics
            U = (k/d).*((pi*L/(4*d))^(-1/3)).*Ra.^(beta); 
            Qs = 2*(S+(1-cc)*Scc)*Kc.*T.*(U./(pi*L*k)).^(1/2); %W
        case 2 % Schubert, using empirical critical Ra (not for onset but for stable convection)
            Qs = 2*(S+(1-cc)*Scc)*Kc.*(T-Ts)/d.*(Ra/Ra_cr).^(1/3); %W
            U = (k/d).*(pi*L/d)*(Ra_cr)^(-2/3).*Ra.^(beta);
        case 3 % Labrosse, 2007: using fixed max plate age & triangle age distribution
            U = L/age_max.*T./T; % not average plate velocity, but subduction rate            
            Qs = (2).*(S)*Kc.*T.*(1./(pi*age_max*k)).^(1/2); %W
        case 4 % Labrosse, 2007: using fixed max plate age 
            U = L/age_max.*(1+cc); % not average plate velocity, but subduction rate 
            %Qs = (2+cc*2/3).*(S+(1-cc)*Scc)*Kc.*T.*(1./(pi*age_max*k)).^(1/2); %W
            Qs = (8/3).*(S)*Kc.*T.*(1./(pi*age_max*k)).^(1/2); %W
            
    end
    
    dl = 2.32*(k.*L./U).^(1/2);% 2.32 from erf

    
%     
%     if Xm==1 % change Xm back to zero as needs be used in calculation for zm and D
%         Xm=0;
%     end
    
    switch etalawin
        case 6 % T already in C for zm
            zm = (z1.*T) + (z2.*Xm.*1e-4) + (z3);
        case 7 % T already in C for zm
            zm = (z1.*T) + (z2.*Xm.*1e-4) + (z3);
        case 99 %          
            zm = (z1.*T) + (z2.*Xm) + (z3);
            zm = max(0,min(zm,300e3));
        otherwise % T needs to be converted from mid-mantle T in K to Tp in C
            Tp_K = T.* exp(-alpha_ra*g*1445e3/Cp); %K potential T calculated from McKenzie and Bickle (1988)
            Tp_C = Tp_K - 273;
            zm = (z1.*(Tp_C)) + (z2.*Xm) + (z3);
            zm = max(0,min(zm,300e3));
    end
    

    % tuned changing R
    R = ((t<=ARt)*ARf1+(t>ARt)*ARf2).*(S+(1-cc)*Scc).*(dl./L).*U.*rhom.*Xp.*Fr ;  %ppm.kg/s
    R_XesAtm = Rxe*((t<=ARt)*ARf1+(t>ARt)*ARf2).*(S+(1-cc)*Scc).*(dl./L).*U.*rhom.*XesAtm.*FrXe ;  %apg. kg/s
    
  
    %changing FD
    
    D0 = (S+(1.-cc)*Scc).*Fd.*(zm./L).*U.*rhom ;  %ppm.kg/s
    D0=(t<=ADt*yr_s)*ADf.*D0+(t>ADt*yr_s).*D0;
    D0_Xe = (S+(1.-cc)*Scc).*Fd.*(zm./L).*U.*rhom;   %apg.kg/s
    D0_Xe=(t<=ADt*yr_s)*ADf.*D0_Xe+(t>ADt*yr_s).*D0_Xe;
    
    D0_p = D0p_pd.*Qc./Qc_pd*FD_p; 
     

    D0_mor = (S+(1-cc)*Scc).*(zm./L).*U.*rhom.*Fd ;
    P0_mor = (S+(1-cc)*Scc).*(zm./L).*U.*rhom;
    

    
    D0_cc = max( ( UrCC.*exp(lamUr.*(t_pd-t)/yr_s)./Umt...
            .*(NG_cc+R_cc) ) - P0_mor,0);


    Hsf_rho_V = Hsf*rho*V;
    
    conH=false;
    
    if conH % if need radiogenic heating to be constant
        H = Hsf_rho_V .* 3.7e-12;
    else     
        % H  = Hsf_rho_V .* sum(H_El,2); %W %2 in sum() sums along each row not column
        H  = rho*V*sum(H_El,2); %W %2 in sum() sums along each row not column
    end
%     

% Fd,Fr,
dynamic_quantity.eta = eta;
dynamic_quantity.Ra = Ra;
dynamic_quantity.U = U;
dynamic_quantity.dl = dl;
dynamic_quantity.zm = zm;
% dynamic_quantity.dM_dt = dM_dt;
dynamic_quantity.R = R;
% dynamic_quantity.D = D;
dynamic_quantity.D0 = D0;
dynamic_quantity.Qs = Qs;
dynamic_quantity.H = H;
%dynamic_quantity.D_Xe = D_Xe;
dynamic_quantity.D0_Xe = D0_Xe;
dynamic_quantity.R_XesAtm = R_XesAtm;
dynamic_quantity.cc = cc;
dynamic_quantity.D0_mor = D0_mor;
dynamic_quantity.D0_cc = D0_cc;
dynamic_quantity.NG_cc = NG_cc;
dynamic_quantity.R_cc = R_cc;
dynamic_quantity.Umt = Umt;
dynamic_quantity.Qc = Qc;
dynamic_quantity.D0_p = D0_p;
