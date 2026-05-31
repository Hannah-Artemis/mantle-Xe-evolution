function [CC_growth_input,Thermal_input,H2O_input,Xe_input,Simulation_input,FIG,FD_MOR,FD_p] = input_degassing(p)



%%
% default input for degassing simulation
% ref original 'multi_runs_FdFr.m'

par_solver_new;
default_input_for_degassing;

 CC_growth_input.time_cc = p.time_cc_input;
        CC_growth_input.kcc = p.kappa_gcc;
        CC_growth_input.krcc = p.kappa_rcc;
        CC_growth_input.Rccs = p.Rs;
        CC_growth_input.Rccp = p.Rp;

    Thermal_input.BSE = p.BSE;
    Thermal_input.para_Scale = [2];
    Thermal_input.Ubse = p.Ubse;
        Thermal_input.Qc_pd = p.Qc_pd_input;
        Thermal_input.dQc = p.dQc_input;
        Thermal_input.etain = p.eta_ref;
        Thermal_input.T_init = p.Ti;

%          switch Thermal_input.BSE
%             case -5 % Javoy-2010-E_BSE (mix of 1/3*EL+2/3*EH)     
%                 [U238bse,U235bse,Thbse,K40bse] = deal(1.26E-08,8.91E-11,3.92E-08,4.30E-08);            
% %             case -4 % Javoy-2010-E_BSE but with lower(normal) K/U
% %                 [U238bse,U235bse,Thbse,K40bse] = deal(1.26E-08,8.91E-11,3.92E-08,1.83E-08);  
%             case -4 % Javoy-2010-E_BSE but with lower(normal) K/U
%                 [U238bse,U235bse,Thbse,K40bse] = deal(1.2641E-08,8.9111E-11,3.9208E-08,1.833E-08);  
%             
%             case -3 % Javoy-1999-E_BSE                  
%                 [U238bse,U235bse,Thbse,K40bse] = deal(1.35E-08,9.54E-11,4.20E-08,4.61E-08);
%             case -2 % Javoy-1999-E_PUM, also used in Labrosse-2007                  
%                 [U238bse,U235bse,Thbse,K40bse] = deal(1.75E-08,1.23E-10,5.54E-08,6.55E-08);
%             case -1 % Jaupart-2015 based on Javoy-2013, considering lower K/U than EN (~1.2e4)                  
%                 [U238bse,U235bse,Thbse,K40bse] = deal(1.79E-08,1.26E-10,6.00E-08,2.58E-08);
%             case 1 % Jellinek-2013_non_chondritic
%                 [U238bse,U235bse,Thbse,K40bse] = deal(1.39E-08,9.80E-11,5.50E-08,1.98E-08);
%             case 2 % LK-2007             
%                 [U238bse,U235bse,Thbse,K40bse] = deal(1.69E-08,1.19E-10,6.30E-08,2.26E-08);
%             case 3 % Javoy-1999-CI_BSE 
%                 [U238bse,U235bse,Thbse,K40bse] = deal(1.95E-08,1.37E-10,6.93E-08,3.21E-08);
%             case 4 % MS-1995             
%                 [U238bse,U235bse,Thbse,K40bse] = deal(1.99E-08,1.40E-10,7.90E-08,2.86E-08);
%             case 5 % Palme and O’Neill (2003)
%                 [U238bse,U235bse,Thbse,K40bse] = deal(2.18E-08,1.54E-10,8.30E-08,3.11E-08);
%          end
        
       
         U238bse = Thermal_input.Ubse*Thermal_input.R0U238_U;
         U235bse = Thermal_input.Ubse*Thermal_input.R0U235_U;
         Thbse =   Thermal_input.Ubse*Thermal_input.R0Th_U;
         K40bse =  Thermal_input.Ubse*Thermal_input.R0K_U*Thermal_input.R0K40_K;

        Thermal_input.h1 = [U238bse;  9.46e-5;    4.4679e9*yr_s     ]; %U238 C30.8 MS13 E1-10.7 E2-6.31
        Thermal_input.h2 = [U235bse;  5.69e-4;    0.704e9*yr_s    ];  %U235 0.22 0.09 0.0757 0.0445
        Thermal_input.h3 = [Thbse;   2.64e-5;    14e9*yr_s       ];    %Th 124 48.6 23.1 10.8
        Thermal_input.h4 = [K40bse;  2.92e-5;    1.25e9*yr_s     ]; %K 36.9 18.7 55.6 36
        Thermal_input.h_El = [Thermal_input.h1,Thermal_input.h2,Thermal_input.h3,Thermal_input.h4];
        
        h_El = Thermal_input.h_El;
        H_El = zeros(4,1);
        H_El(1) = (h_El(1,1)-U238CC*mcc_val/rho/V).*h_El(2,1);
        H_El(2) = (h_El(1,2)-U235CC*mcc_val/rho/V).*h_El(2,2);
        H_El(3) = (h_El(1,3)-ThCC*mcc_val/rho/V).*h_El(2,3);
        H_El(4) = (h_El(1,4)-KCC*mcc_val/rho/V).*h_El(2,4);
        H_pd = sum(H_El)*rho*V;
        Thermal_input.Hsf = 1+Thermal_input.Qc_pd/H_pd;
       

    H2O_input.Frw = p.Fr_w;
        H2O_input.Xm_init = p.Xm_init;
        H2O_input.Ms_init = p.Ms_init;

    Xe_input.SwXe = p.SwXe;
        Xe_input.Xe_init = p.Xe_init;
        Xe_input.FrXe = p.Fr_Xe;

        if Thermal_input.BSE > 0
            Xer128=Xercc; Xer131=Xe131rcc; Xer132=Xe132rcc; Xer134=Xe134rcc; Xer136=Xe136rcc;
        else
            Xer128=Xeren; Xer131=Xe131ren; Xer132=Xe132ren; Xer134=Xe134ren; Xer136=Xe136ren;
        end
        Xe_input.XeAtm_init    = Xer128*Xe_input.Xe_init ; % atmospheric Xe in mantle ->initial 128Xe in mantle 3.5054e-9
        Xe_input.Xe136_init =  Xer136*Xe_input.Xe_init ;%Xe136rcc*130Xe_ini 1.3737e-8
        Xe_input.Xe131_init = Xer131*Xe_input.Xe_init;%Xe136rcc*130Xe_ini 
        Xe_input.Xe132_init =  Xer132*Xe_input.Xe_init ;%Xe132rcc*130Xe_ini
        Xe_input.Xe134_init = Xer134*Xe_input.Xe_init ;%Xe136rcc*130Xe_ini
        
        Xe_input.lamUr=log(2)/(Thermal_input.h1(3))*yr_s;
        Xe_input.Ur_init=U238bse*exp(Xe_input.lamUr*t_pd/yr_s)/238*NA;%atoms/g
        Xe_input.lamPu=8.6643e-9;%yr-1
        Xe_input.R0Pu_Ur=0.0068;%initial Pu/U 0.0068 (atomic ratio)
        Xe_input.Pu_init=Xe_input.Ur_init*Xe_input.R0Pu_Ur;%atoms/g

    FD_MOR = p.Fd_mor; 
    FD_p = p.Fd_p;



