

%%
% default input for degassing simulation
% ref original 'multi_runs_FdFr.m'

par_solver_new;


CC_growth_input.time_cc = 0.3e9*yr_s;% start point for CC growth
    CC_growth_input.kcc= 4.0756;
    %finished 1Ga before:1.5322; 2Ga before:2.2992; 3Ga before 4.0756; 
    % 2.5Ga:2.9432
    CC_growth_input.Rccs=10e22;%0-1e23 kg/Ga
    CC_growth_input.Rccp=1e22;%0-2e22 kg/Ga
    CC_growth_input.krcc=1.0;%(-3,3)



Thermal_input.T_init        = [ 2500    ]; % Initial Average Mantle Temperature (K) 1883K for Tp with Cp=1250
    Thermal_input.BSE        = [ -4    ]; 
    Thermal_input.Ubse        = [ 15e-9    ]; 
    Thermal_input.Qc_pd        = [ 15e12    ]; %W
    Thermal_input.dQc        = [ 2e12    ]; %W
    % Thermal_input.Hsf      = [ 1  ]; % Scaling for Radiogenic Heating 1.706;2.1(for UMS model);2.8(for UEH3)
    Thermal_input.para_Scale = [2];% 1 Crowley; 2 Schubert; 3 Labrosse 5 Schubert but with Tp
    Thermal_input.etain      = 7.5e14; % 10=Ze14 12=Ze15 14=Ze16 16=Ze17 18=Ze18 20=Ze19 % Z=8.73 for eta_pd=1e22; Z=3.77 for eta_pd=5e21
    Thermal_input.etalawin = 100;
    Thermal_input.beta     = [ 2/3    ];
    % Prescribed Plate Velocity (i.e. NOT Ra based calculations) 
    Thermal_input.pU = logical(0); 
        Thermal_input.U_cm_yr = [ 0 ] ; %cm/yr
    
%     switch Thermal_input.BSE
%             case -5 % Javoy-2010-E_BSE (mix of 1/3*EL+2/3*EH)     
%                 [U238bse,U235bse,Thbse,K40bse] = deal(1.26E-08,8.91E-11,3.92E-08,4.30E-08);            
%             case -4 % Javoy-2010-E_BSE but with lower(normal) K/U
%                 [U238bse,U235bse,Thbse,K40bse] = deal(1.26E-08,8.91E-11,3.92E-08,1.83E-08);  
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
%     end

    Thermal_input.R0U238_U = 0.9927;
        Thermal_input.R0U235_U = 1-Thermal_input.R0U238_U;
        Thermal_input.R0Th_U = 4;% 4[GK-CI]; 3.08; 3.0-3.8[2014-Javoy-EN]
        Thermal_input.R0K_U = 1.27e4;% 1.27e4[GK];1.21e4
        Thermal_input.R0K40_K = 1.28e-4;% 1.28e-4[GK];1.19e-4 
       
         U238bse = Thermal_input.Ubse*Thermal_input.R0U238_U;
         U235bse = Thermal_input.Ubse*Thermal_input.R0U235_U;
         Thbse =   Thermal_input.Ubse*Thermal_input.R0Th_U;
         K40bse =  Thermal_input.Ubse*Thermal_input.R0K_U*Thermal_input.R0K40_K;

    Thermal_input.h1 = [U238bse;  9.46e-5;    4.4679e9*yr_s     ]; %U238 C30.8 MS13 E1-10.7 E2-6.31
    Thermal_input.h2 = [U235bse;  5.69e-4;    0.704e9*yr_s    ];  %U235 0.22 0.09 0.0757 0.0445
    Thermal_input.h3 = [Thbse;   2.64e-5;    14e9*yr_s       ];    %Th 124 48.6 23.1 10.8
    Thermal_input.h4 = [K40bse;  2.92e-5;    1.25e9*yr_s     ]; %K 36.9 18.7 55.6 36
    Thermal_input.h_El = [Thermal_input.h1,Thermal_input.h2,Thermal_input.h3,Thermal_input.h4];

    H_El = zeros(4,1);h_El = Thermal_input.h_El;
        H_El(1) = (h_El(1,1)-U238CC*mcc_val/rho/V).*h_El(2,1);
        H_El(2) = (h_El(1,2)-U235CC*mcc_val/rho/V).*h_El(2,2);
        H_El(3) = (h_El(1,3)-ThCC*mcc_val/rho/V).*h_El(2,3);
        H_El(4) = (h_El(1,4)-KCC*mcc_val/rho/V).*h_El(2,4);
        H_pd = sum(H_El)*rho*V;
    Thermal_input.Hsf = 1+Thermal_input.Qc_pd/H_pd;


H2O_input.Xm_init       = [ 320 ]; % Initial Mantle Water Concentration (ppm) 2500    2500    2500    2500    2500
    H2O_input.Ms_init       = [ 1.5 ] ; % Initial Surface Ocean Mass (multiples of present day Earth ocean mass, M_oc)
    H2O_input.Frw       = 0.25;
    H2O_input.Fd_sw        = [ logical(0)  ]; %R=0, D remains when Ms=0; % degassing efficiency switches to 1 after 3 OM in the mantle
        H2O_input.Fd_new   = [      1      ];
    H2O_input.f_WI        = [ 302  ]; % ppm


%IN.Fd       =  Fd;%[ 0 0.001 0.01 0.1 1 ]; %[ 0:0.001:0.009,0.01:0.01:0.09,0.1:0.1:1 ]; %[ 0 0.001 0.01 0.1 1 ];
    


% 130Xe Tracking
Xe_input.SwXe         = [ logical(1) ];
    Xe_input.Xe_init       = [ 3.2e8   ]; % primordial mantle 130Xe ppm 6.91e-9(=3.2e7atom/g)
    Xe_input.XesAtm_init   = [ 8e6   ]; % Xe in regassing (in hydrous mineral) ppm 1.262e-9 (=5.85e6atoms/g)
    Xe_input.FrXe = 0.2;

    Xe_input.Xes_init      = [ Xerpa*(1+atmf*(130-128))    ]; % primordial surface Xe -> initial ATM 128Xe/130Xe 0.5083
    Xe_input.Xe136r130atm_init = [ Xe136rpa*(1+atmf*(130-136)) ]; %Xe136rpa after back to initial 1.6668 
    Xe_input.Xe131r130atm_init = [ Xe131rpa*(1+atmf*(130-131)) ]; %Xe136rpa after back to initial
    Xe_input.Xe132r130atm_init = [ Xe132rpa*(1+atmf*(130-132)) ]; %Xe136rpa after back to initial
    Xe_input.Xe134r130atm_init = [ Xe134rpa*(1+atmf*(130-134)) ]; %Xe136rpa after back to initial

    
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

    %Xe_input.Xe128    = [ 3.5054e-9  ]; % primordial mantle 128Xe ppm   
    %Xp:XesAtm in mass for serp = 1.63E+12
    Xe_input.lamUr=log(2)/(Thermal_input.h1(3))*yr_s;
    Xe_input.Ur_init=U238bse*exp(Xe_input.lamUr*t_pd/yr_s)/238*NA;%atoms/g
    Xe_input.lamPu=8.6643e-9;%yr-1
    Xe_input.R0Pu_Ur=0.0068;%initial Pu/U 0.0068 (atomic ratio)
    Xe_input.Pu_init=Xe_input.Ur_init*Xe_input.R0Pu_Ur;%atoms/g

    

% Solver Settings
Simulation_input.end      = 1; % run solver for a number of earth lifetimes i.e. 1 => 4.6Gyrs 2 => 2x4.6Gyrs etc...
    Simulation_input.etalaw   = Thermal_input.etalawin; % Viscosity Law (see eta_calc.m for more details)
    % Simulation_input.dynsol   = [ 4      4       4    ]; % Solver for Dynamic Quantities (unless investigati
    % Simulation_input.RK       = [ 0      0       0    ]; % ODE solver: 1=Runge-Kutta (with time lag) 0=ode15s 
    
%     Simulation_input.nsteps  = Simulation_input.end*[ 46001 ]; % number of timesteps
%     t_fin = Simulation_input.end*t_pd;
%     t = linspace(0,t_fin,Simulation_input.nsteps);
%     Simulation_input.time_series = t';
    t_set = [0:0.0001:2,2.001:0.001:2.6,2.601:0.001:4.6]/4.6*t_pd;
    Simulation_input.nsteps=length(t_set);
    Simulation_input.time_series = t_set';
    %Simulation_input.cases   = length(IN.Fd); % number of cases to compare (no. of columns above)
    % Mixing Parameters
    Simulation_input.mixing       = [ logical(0) ];
        Simulation_input.mixlag   = [ 2e9        ]; % lag time for mixing in years 
        Simulation_input.Rapd     = [ 1e6        ]; %present day Rayleigh number for mixing lag calcaultion
        Simulation_input.varmix   = [ 0          ]; % 1 = variable mixing lag scaled with Ra, 0 = constant mixing lag
        Simulation_input.m        = [ 2/3        ];
    
    % Post Processing (for mixing only)
    Simulation_input.postpro = 0;



% Figure Settings

FIG.make        = [1]; % generate figures; see fig_plot for details 1=subplots 6=individual figures (for posters etc.)
    FIG.makepp = [10,11];
    FIG.save        = [1]; % save figure to FIG.savename (only for indivdual figures, subplots saved as default)
    FIG.savename    = '0.5OMXeall_Rlimited_apg/newtest'; % '/Volumes/KCHARDDRIVE/PhD/#ParamModel_170518/#mixingpaper/#analysis_new_long/#10Gyrs/O0_sca/O0sca';data saved in working directory as FIG.savename.mat file
    FIG.fntsz       = 20; % font size (suggested 36 for individual figures, 20 for subplots)
    FIG.comp        = [0]; % 0=no plot comparison 1=Crowley2011 2=San2011 3='Corrected' Crowley
    FIG.lim         = [0 Simulation_input.end*4.6]; % time axis range (Gyrs)
    FIG.calib       = [1]; % show present day calibration guides = 1 
    FIG.lc          = [ 1 0 0; 0 1 1; 0 0 1; 1 0 1; 0 1 0; 1 1 0;...
                        1 0 0; 0 1 1; 0 0 1; 1 0 1; 0 1 0; 1 1 0;...
                        1 0 0; 0 1 1; 0 0 1; 1 0 1; 0 1 0; 1 1 0;...
                        1 0 0; 0 1 1; 0 0 1; 1 0 1; 0 1 0; 1 1 0;...
                        1 0 0; 0 1 1; 0 0 1; 1 0 1; 0 1 0; 1 1 0];
                        %[ 192 192 192; 255 128 0; 0 204 0;192 192 192; 255 128 0; 0 204 0 ]./255; 
                        %[ 1 0 0; 0 1 1; 0 0 1; 1 0 0; 0 1 1; 0 0 1;1 0 1; 0 1 0; 1 1 0]; % matrix of rgb colours for each case
    FIG.ls          = char('-','-','-','-','-','-',...
                            '-','-','-','-','-','-',...
                            '-','-','-','-','-','-',...
                            '-','-','-','-','-','-',...
                            '-','-','-','-','-','-');
                        %char('-','--',':','-','--',':'); 
                        %char('-','-','-','--','--','--'); % line styles for each case



