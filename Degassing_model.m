function output = Degassing_model(Simulation_input,Thermal_input,CC_growth_input,H2O_input,Xe_input,FD_MOR,FD_p)


fin = Simulation_input.end;

% thermal model related parameters
T_init = Thermal_input.T_init;
BSE = Thermal_input.BSE;
Ubse = Thermal_input.Ubse;
h_El = Thermal_input.h_El;
Qc_pd = Thermal_input.Qc_pd;
dQc =Thermal_input.dQc;
Hsf = Thermal_input.Hsf;
etain = Thermal_input.etain;
etalawin =  Thermal_input.etalawin;
para_Scale = Thermal_input.para_Scale;
beta = Thermal_input.beta;
% pU = Thermal_input.pU; % prescribed velocity
U_cm_yr = Thermal_input.U_cm_yr;



% CC grwoth model related parametes
time_cc = CC_growth_input.time_cc;
kcc = CC_growth_input.kcc;
krcc = CC_growth_input.krcc;
Rccs = CC_growth_input.Rccs;
Rccp = CC_growth_input.Rccp;


%
SwXe = Xe_input.SwXe;
Xe_init = Xe_input.Xe_init;
Xes_init = Xe_input.Xes_init;
XeAtm_init = Xe_input.XeAtm_init;
XesAtm_init = Xe_input.XesAtm_init;
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
R0Pu_Ur=Xe_input.R0Pu_Ur;
lamUr = Xe_input.lamUr;
lamPu = Xe_input.lamPu;
FrXe = Xe_input.FrXe;


Fd=FD_MOR;




%% Input Parameters
% tic

par_solver_new;
%t_fin = fin*4.6e9*yr_s;
t_fin = fin*t_pd;
% t = linspace(0,t_fin,nsteps); % define time series
% t=t';
t=Simulation_input.time_series;
% nsteps = Simulation_input.nsteps;
nsteps = length(t);
deltat = t(2)-t(1);

%prescribed_U = pU;
U = U_cm_yr / 100 / yr_s; % cm/yr -> m/s



%% Solving


% Euler method
  

%initialize
Xe = [Xe_init;zeros(nsteps-1,1)];
Xes = [Xes_init;zeros(nsteps-1,1)]; % 128Xe/130Xe in atmosphere
XeAtm = [XeAtm_init;zeros(nsteps-1,1)]; % 128Xe in mantle
XesAtm = [Xes_init;zeros(nsteps-1,1)]; % 128Xe/130Xe in ?
Xe136r130atm = [Xe136r130atm_init;zeros(nsteps-1,1)];
Xe136 = [Xe136_init;zeros(nsteps-1,1)];
Xe131r130atm = [Xe131r130atm_init;zeros(nsteps-1,1)];
Xe131 = [Xe131_init;zeros(nsteps-1,1)];
Xe132r130atm = [Xe132r130atm_init;zeros(nsteps-1,1)];
Xe132 = [Xe132_init;zeros(nsteps-1,1)];
Xe134r130atm = [Xe134r130atm_init;zeros(nsteps-1,1)];
Xe134 = [Xe134_init;zeros(nsteps-1,1)];
Put= [Pu_init;zeros(nsteps-1,1)];
Urt= [Ur_init;zeros(nsteps-1,1)];          

% derive temporal CC & T evolution
[Tp,Tm,Qs,Qc,H,Ra,eta] = ...
Thermal_model(t,T_init,etain,Qc_pd,dQc,...
            BSE,Ubse,kcc,time_cc,para_Scale,Hsf);
T=Tm;

% dynamic quantities
dynamic_quantity = ...
        dynamics_solver_new(t,0,Tm,Thermal_input,CC_growth_input,H2O_input,Xe_input,FD_MOR,FD_p);
% end

eta=dynamic_quantity.eta;
Ra=dynamic_quantity.Ra;
U=dynamic_quantity.U; % velocity
dl=dynamic_quantity.dl;
zm=dynamic_quantity.zm;
%dM_dt=dynamic_quantity.dM_dt;
R=dynamic_quantity.R;
D0=dynamic_quantity.D0;
Qs=dynamic_quantity.Qs;
H=dynamic_quantity.H;
D0_Xe=dynamic_quantity.D0_Xe;
R_XesAtm=dynamic_quantity.R_XesAtm;
cc=dynamic_quantity.cc;
D0_mor=dynamic_quantity.D0_mor;
D0_cc=dynamic_quantity.D0_cc;
NG_cc=dynamic_quantity.NG_cc;
R_cc=dynamic_quantity.R_cc;
Umt=dynamic_quantity.Umt;
Qc=dynamic_quantity.Qc;
D0_p=dynamic_quantity.D0_p;


Urt=Umt;% in atoms/g

for i = 2:nsteps

    % SPu=-lamPu/yr_s*Put(i-1);
    SPu=0;
    kPu=lamPu/yr_s...
        +UrCC*exp(lamUr*(t_pd-t(i-1))/yr_s)/Urt(i-1)*NG_cc(i-1)/V/rho;
    ekPu=exp(-kPu*(t(i)-t(i-1)));

   % Put(i) = Put(i-1) + (t(i)-t(i-1)) * (SPu-Put(i-1)*kPu);
   % Put(i) = ( Put(i-1) + (t(i)-t(i-1))*SPu ) ./ ( 1 + (t(i)-t(i-1))*kPu );
   Put(i) = Put(i-1).*ekPu + (SPu./max(k,eps)) .* (1 - ekPu);

    % 130Xe Dxecc=max(Dxecc0-Dxe,0)
    Xe(i) = Xe(i-1) + (t(i)-t(i-1)) * ( (R_XesAtm(i-1)-D0_mor(i-1)*Xe(i-1)) / (rho*V) - DXecc*Xe(i-1)*D0_cc(i-1)/V/rho - DXep*Xe(i-1)*D0_p(i-1)/V/rho ); % actually dx/dt
    XesAtm = 0; % XesAtm
   
    % ATM iXe/130Xe
    if t(i)<2.6*1e9*yr_s
        Xes(i) = Xes(i-1) + (t(i)-t(i-1)) * (-Xerpa*(130-128)*39/1000*1/(2.6*1e9*yr_s)); % 128Xe/130Xe 
        Xe131r130atm(i) = Xe131r130atm(i-1) + (t(i)-t(i-1)) * (-Xe131rpa*(130-131)*39/1000*1/(2.6*1e9*yr_s));
        Xe132r130atm(i) = Xe132r130atm(i-1) + (t(i)-t(i-1)) * (-Xe132rpa*(130-132)*39/1000*1/(2.6*1e9*yr_s));
        Xe134r130atm(i) = Xe134r130atm(i-1) + (t(i)-t(i-1)) * (-Xe134rpa*(130-134)*39/1000*1/(2.6*1e9*yr_s));
        Xe136r130atm(i) = Xe136r130atm(i-1) + (t(i)-t(i-1)) * (-Xe136rpa*(130-136)*39/1000*1/(2.6*1e9*yr_s));
    else
        Xes(i) = Xes(i-1) + (t(i)-t(i-1)) * 0;
        Xe131r130atm(i) = Xe131r130atm(i-1);
        Xe132r130atm(i) = Xe132r130atm(i-1);
        Xe134r130atm(i) = Xe134r130atm(i-1);
        Xe136r130atm(i) = Xe136r130atm(i-1);
    end    

  

    % 128Xe in the mantle
    XeAtm(i) = XeAtm(i-1) + (t(i)-t(i-1)) * ( (R_XesAtm(i-1)*Xes(i-1)...
        -D0_mor(i-1)*XeAtm(i-1)) / (rho*V) ...
        - DXecc*XeAtm(i-1)*D0_cc(i-1)/V/rho - DXep*XeAtm(i-1)*D0_p(i-1)/V/rho ); 

    par_timestep=(t(i)-t(i-1));

    % 131Xe in mantle
    Xe131(i) = Xe131(i-1) + (t(i)-t(i-1)) * ( ( R_XesAtm(i-1)*Xe131r130atm(i-1)- ...
        ( D0_mor(i-1) + DXecc*D0_cc(i-1) + DXep*D0_p(i-1) )...
        *(Xe131(i-1)+(lamPu/yr_s*YPu131*Put(i-1)+lamUr/yr_s*YUr131*Urt(i-1))*par_timestep)) ...
        / (rho*V) ...
        +lamPu/yr_s*YPu131*Put(i-1)+lamUr/yr_s*YUr131*Urt(i-1) ); 

 
    % 132Xe in mantle
    Xe132(i) = Xe132(i-1) + (t(i)-t(i-1)) * ( ( R_XesAtm(i-1)*Xe132r130atm(i-1) - ...
        ( D0_mor(i-1) + DXecc*D0_cc(i-1) + DXep*D0_p(i-1))...
        *(Xe132(i-1)+(lamPu/yr_s*YPu132*Put(i-1)+lamUr/yr_s*YUr132*Urt(i-1))*par_timestep) ) ...
        / (rho*V) ...
        +lamPu/yr_s*YPu132*Put(i-1)+lamUr/yr_s*YUr132*Urt(i-1) ); 

  
    % 134Xe in mantle
    Xe134(i) = Xe134(i-1) + (t(i)-t(i-1)) * ( ( R_XesAtm(i-1)*Xe134r130atm(i-1)-...
         ( D0_mor(i-1) + DXecc*D0_cc(i-1) + DXep*D0_p(i-1))...
        *( Xe134(i-1)+(lamPu/yr_s*YPu134*Put(i-1)+lamUr/yr_s*YUr134*Urt(i-1))*par_timestep) ) ...
        / (rho*V)...
        +lamPu/yr_s*YPu134*Put(i-1)+lamUr/yr_s*YUr134*Urt(i-1) ); 

    % 136Xe in mantle     
    Xe136(i) = Xe136(i-1) + (t(i)-t(i-1)) * ( ( R_XesAtm(i-1)*Xe136r130atm(i-1) - ...
           ( D0_mor(i-1) + DXecc*D0_cc(i-1) + DXep*D0_p(i-1) ) ...
          *( Xe136(i-1)+(lamPu/yr_s*YPu*Put(i-1)+lamUr/yr_s*YUr*Urt(i-1))*par_timestep) ) ...
          / (rho*V)...
          +lamPu/yr_s*YPu*Put(i-1)+lamUr/yr_s*YUr*Urt(i-1) ); 


end




%% Save to Output
%output = zeros(1,37);

output.t_Ga = t/(1e9*yr_s);   %Ga
output.T = T;                 %K Tavg
output.U = U*100*yr_s;       %cm/yr
output.eta = eta;            %Pas
output.Ra = Ra; 
output.dl = dl/1e3;          %km

output.zm = zm/1e3;          %km
output.Xp = Xp;              %ppm
output.Fd = Fd;              %dimensionless
output.Fd_mor = FD_MOR;
output.Fd_plume = FD_p;
output.R = R*1e-6*yr_s;      %kg/yr
output.D0 = D0*1e-6*yr_s;      %kg/yr
output.Qs = Qs;              %W
output.H = H;                  %W
output.Qc = Qc;
output.Hsf = Hsf;      
output.etalaw = etalawin; 
output.etain = etain;
output.Ur = H./Qs;           % Urey Ratio
output.Nu = (Qs.*d)./(Kc.*T.*S); % Nusselt Number
output.timesteps  = nsteps; % number of timesteps
output.Qc = Qc;
output.D0_p = D0_p;
output.cc = cc;
output.D0_mor = D0_mor; % mor degassing kg/s
output.Pmor = D0_mor*yr_s/FD_MOR; % mor mantle processing rate kg/yr
output.D0_cc = D0_cc; % cc generation degassing kg/s

output.SwXe = SwXe;
%if SwXe
    output.Fr_Xe = FrXe;
    output.Xei = Xe_init;
    output.Xed = XesAtm_init;
    output.Urini = Ur_init;
    output.R0Pu_Ur = R0Pu_Ur;
    output.Xe = Xe; % xenon concentration (ppm)
    output.D0_Xe = D0_Xe; % xenon degassing atoms/g*kg/s    
    output.NG_cc = NG_cc; % cc net growth kg/s
    output.Umt = Umt;
    output.R_cc = R_cc; % cc recycling kg/s
    output.Xes = Xes; % surface Xe concentration (ppm)    
    output.R_Xe = R_XesAtm; % atmospheric xenon regassing
    output.Xe_Atm = XeAtm; % atmospheric mantle Xe concentration (ppm)
    output.Xes_Atm = XesAtm; % atmospheric surface Xe concentration (ppm)
    output.Xe136r130atm = Xe136r130atm;
    output.Xe136 = Xe136;
    output.Xe131r130atm = Xe131r130atm;
    output.Xe131 = Xe131;
    output.Xe132r130atm = Xe132r130atm;
    output.Xe132 = Xe132;
    output.Xe134r130atm = Xe134r130atm;
    output.Xe134 = Xe134;
    output.Put=Put;
    output.Urt=Urt;

end