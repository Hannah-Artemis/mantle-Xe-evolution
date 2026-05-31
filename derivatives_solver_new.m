function [dX] = derivatives_solver_new(t,X,Xm_lag,Simulation_input,Thermal_input,CC_growth_input,H2O_input,Xe_input,FD_MOR,FD_p)

% X(1) = Mantle Water Content
% X(2) = Average Mantle Temperature
% X(3) = Surface Ocean Mass
% 
% global etalawin
% global para_Scale
% % global dynsol
% global etain
% global Fr
% global Fd
% global Hsf
% global prescribed_U
% global U   
% global Fd_sw
% global Fd_new
% global beta
% global Ms_init
% global Fd_sw
% global XePu %136Xe from Pu fission
% global XeUr %136Xe from U fission


% CC_growth_input=params.CC_growth_input;
% Thermal_input=params.Thermal_input;
% H2O_input=params.H2O_input;
% Xe_input=params.Xe_input;

par_solver_new;

nsteps = Simulation_input.nsteps;
fin = Simulation_input.end;
mixing = Simulation_input.mixing;
mlag_yrs = Simulation_input.mixlag;
varmix = Simulation_input.varmix;
m = Simulation_input.m;
Rapd = Simulation_input.Rapd;
par_timestep = t_pd./(Simulation_input.nsteps./Simulation_input.end);


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
prescribed_U = Thermal_input.pU;
U_cm_yr = Thermal_input.U_cm_yr;



% CC grwoth model related parametes
time_cc = CC_growth_input.time_cc;
kcc = CC_growth_input.kcc;
krcc = CC_growth_input.krcc;
Rccs = CC_growth_input.Rccs;
Rccp = CC_growth_input.Rccp;

%
Xm_init = H2O_input.Xm_init;
Ms_init = H2O_input.Ms_init;
Fr = H2O_input.Frw;
Fd_sw = H2O_input.Fd_sw;
Fd_new = H2O_input.Fd_new;
f_WI = H2O_input.f_WI;

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
R0Pu_Ur = Xe_input.R0Pu_Ur;
lamUr = Xe_input.lamUr;
lamPu = Xe_input.lamPu;
FrXe = Xe_input.FrXe;


Fd=FD_MOR;
%Fd_p =



%disp(length(t));
% if number or arguments in is greater than or equal to 3 and 
if nargin >= 3 && ~isempty(Xm_lag)>0
    Xm = Xm_lag;
else
    Xm  = X(1); % Xm_lag;;
end

% Xm  = Xm_lag;
T   = X(2);
Ms  = X(3);
% Tc  = X(4);

if (length(X)>3)
    Xe      = X(4); % primordial Xe in mantle
    Xes     = X(5); % primordial Xe at surface
    XeAtm   = X(6); % atmospheric Xe in mantle
    XesAtm  = X(7); % atmospheric Xe at surface
    Xe136r130atm = X(8); % atmospheric 136Xe/130Xe
    Xe136 = X(9);
    Xe131r130atm = X(10); % atmospheric 136Xe/130Xe
    Xe131 = X(11);
    Xe132r130atm = X(12); % atmospheric 136Xe/130Xe
    Xe132 = X(13);
    Xe134r130atm = X(14); % atmospheric 136Xe/130Xe
    Xe134 = X(15);
    Put=X(16);
    Urt=X(17);

else
    Xe = 0;
    XesAtm = 0;
end

% switch dynsol
%     case 1
%         [eta,Ra,U,dl,zm,dM_dt,Xp,Fd,R,D,Qs,H,Fc,delb,Qcmb] = dynamics_solver(t,Xm,T,Tc,etalawin,etain,Fr,Fd);
%     case 2
%         [eta,Ra,U,dl,zm,dM_dt,Xp,Fd,R,D,Qs,H] = dynamics_solver_paper_results(t,Xm,T,etalawin,etain);
%     case 3
%         [eta,Ra,U,dl,zm,dM_dt,Xp,Fd,R,D,Qs,H] = dynamics_solver_NuRa(t,Xm,T,etalawin);
%     case 4

        if prescribed_U
            
            [eta,Ra,U,dl,zm,dM_dt,Fd,R,D,Qs,H,Fr,D_Xe,R_XesAtm] = ...
                dynamics_solver_U(t,Xm,T,etalawin,etain,Fr,Fd,Xe,Hsf,XesAtm,U);
         else
                
%                 [eta,Ra,U,dl,zm,dM_dt,Fd,R,D,Qs,H,Fr,D_Xe,R_XesAtm,cc,D0_mor,D0_cc,NG_cc,R_cc,Umt,Qc,D0_p] = ...
%             dynamics_solver_paper_maths(t,Xm,T,etalawin,etain,Fr,Fd,Xe,Hsf,XesAtm,Fd_sw,Fd_new,beta,para_Scale);

%              dynamic_quantity = ...
%                     dynamics_solver_new(t,Xm,T,Xe,Thermal_input,CC_growth_input,H2O_input,Xe_input,FD_MOR,FD_p);
              dynamic_quantity = ...
                    dynamics_solver_all(t,Xm,T,Xe,Thermal_input,CC_growth_input,H2O_input,Xe_input,FD_MOR,FD_p);


        end
        
        
%     otherwise 
%         error(['Unknown dynsol: ' dynsol])  
% end




eta=dynamic_quantity.eta;
Ra=dynamic_quantity.Ra;
U=dynamic_quantity.U;
dl=dynamic_quantity.dl;
zm=dynamic_quantity.zm;
dM_dt=dynamic_quantity.dM_dt;
R=dynamic_quantity.R;
%D=dynamic_quantity.D;
D0=dynamic_quantity.D0;
Qs=dynamic_quantity.Qs;
H=dynamic_quantity.H;
%D_Xe=dynamic_quantity.D_Xe;
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

% % from Honing & Spohn (2016)
% Vcore = (4*pi/3)*(r_c^3); %1.75e20; %m3
% Ac = 4.*pi.*r_c.^2;%1.514e14; %m.sq 
% cpcrho = 3.6e6; %J.K-1.m-3 


% %%
% %Pu & U
%  persistent tm;
%  persistent cc_test;
%  persistent deltat;
% persistent Pu;
% persistent Pum;
% persistent Ur;
% persistent Urm;
% persistent kk;
%  if isempty(tm)
%      tm=0;
%      cc_test=0;
%      deltat=0;%kk=1;
%  else     
%      deltat=t-tm;
%      tm=t;
%      cc_test=cc_test+deltat*(1/1e9/yr_s*((t<300e6*yr_s)*0+(t>=300e6*yr_s)...
%         .*(2*kcc*exp(kcc*(t-300e6*yr_s)/yr_s/1e9)...
%         .*(exp(kcc*(t_pd-time_cc))+1)...
%         ./((exp(kcc*(t_pd-time_cc))-1)*((exp(-kcc*(300e6*yr_s-t)/yr_s/1e9)+1)).^2))));
%      %kk=kk+1;
%  end
%  scatter(t/yr_s/1e9,cc_test,'k');hold on;
% %disp(deltat);
% if isempty(Ur)
%     Urm=Urini;Ur=Urm;
% else
%     Ur=Urm*exp(-lamUr*abs(deltat)/yr_s)-UrCC*((t<300e6*yr_s)*0+(t>300e6*yr_s)*(abs(deltat)/4.3e9/yr_s))*mcc/V/rho;
%     Urm=Ur;
% end
% if isempty(Pu)
%     Pum=Puini;Pu=Pum;
% else
%     Pu=Pum*exp(-lamPu*abs(deltat)/yr_s)-UrCC*Pum/Urm*((t<300e6*yr_s)*0+(t>300e6*yr_s)*(abs(deltat)/4.3e9/yr_s))*mcc/V/rho;
%     Pum=Pu;
% end
% ttt(kk)=t;
% nUri(kk)=Ur;
% nPu(kk)=Pu;
% %disp(t);disp(Ur)
% %scatter(t/yr_s/1e9,deltat/yr_s/1e6,'k');hold on;
% %scatter(t/yr_s/1e9,Pu,'k');hold on;
% %disp(outputUri(kk));
% 
% XePu=Pum*(1-exp(-lamPu*abs(deltat)/yr_s))*YPu;
% XeUr=Urm*(1-exp(-lamUr*abs(deltat)/yr_s))*YUr;
% %scatter(t/yr_s/1e9,(XePu+XeUr)/Xe136,'k');hold on;
% scatter(t/yr_s/1e9,Put,'k');hold on;

%%
%dX/dt
% dX(1) = (R-D)/(rho*V)-Dwcc*max( (UrCC*exp(lamUr*(t_pd-t)/yr_s)*X(1)/Urt...
%         *(1/1e9/yr_s*( (t<time_cc)*0+(t>=time_cc)... %net CC growth; Urcc and Urt in atoms/g; X(1) in ppm; thus Xmdot = ppm/s
%         .*( ( 2*kcc*exp(kcc*(t-time_cc)/yr_s/1e9)...
%         .*(exp(kcc*(t_pd-time_cc)/yr_s/1e9)+1)...
%         ./((exp(kcc*(t_pd-time_cc)/yr_s/1e9)-1)*((exp(-kcc*(time_cc-t)/yr_s/1e9)+1)).^2) ) ...
%         +( Rccs/mcc+ (Rccp-Rccs)/mcc.*(1-exp(-krcc*(t-time_cc)/yr_s/1e9))./(1-exp(-krcc*(t_pd-time_cc)/yr_s/1e9)) )  ) ) )... %CC recycling
%         *mcc)-D,0)/V/rho; % %Xmdot = ppm/s 

%dX(1) = (R-D)/(rho*V)-Dwcc*X(1)*D0_cc/V/rho-Dwp*D0_p*X(1)/V/rho; % %Xmdot = ppm/s 
dX(1) = (R-D0_mor*X(1))/(rho*V)-Dwcc*X(1)*D0_cc/V/rho-Dwp*D0_p*X(1)/V/rho; 

%dX(2) = (-Qs+H+Qc)/(Cp*rho*V); %Tdot = K/s %Qcmb
dX(2) = (-Qs+H+Qc)/(Cp*rho*V).*exp(alpha_ra*g*1445e3/Cp); %Tdot = K/s %Qcmb
% dX(3) = dM_dt*(1/(M_oc*1e6))*((Xm*Fd)-(Xp*Fr)); %dMs/dt = (Mw/Moc)*(1/s)
%dX(3) =  (1/(M_oc*1e6))*(D-R); %1e-6*rho*V/M_oc = 0.0023
% dX(3) =  (1/(M_oc*1e6))*(D-R+Dwcc*max( (UrCC*exp(lamUr*(t_pd-t)/yr_s)*X(1)/Urt...
%         *(1/1e9/yr_s*( (t<time_cc)*0+(t>=time_cc)... %net CC growth; Urcc and Urt in atoms/g; X(1) in ppm; thus Xmdot = ppm/s
%         .*( ( 2*kcc*exp(kcc*(t-time_cc)/yr_s/1e9)...
%         .*(exp(kcc*(t_pd-time_cc)/yr_s/1e9)+1)...
%         ./((exp(kcc*(t_pd-time_cc)/yr_s/1e9)-1)*((exp(-kcc*(time_cc-t)/yr_s/1e9)+1)).^2) ) ...
%         +( Rccs/mcc+ (Rccp-Rccs)/mcc.*(1-exp(-krcc*(t-time_cc)/yr_s/1e9))./(1-exp(-krcc*(t_pd-time_cc)/yr_s/1e9)) )  ) ) )... %CC recycling
%         *mcc)-D,0)); 
%dX(3) =  (1/(M_oc*1e6))*(D-R+Dwcc*X(1)*D0_cc+Dwp*X(1)*D0_p); 
dX(3) =  (1/(M_oc*1e6))*(D0_mor*X(1)-R+Dwcc*X(1)*D0_cc+Dwp*X(1)*D0_p); 
% dX(4) = (-Ac*Fc)/(cpcrho*Vcore);



if (length(X)>3)

%     % Primodial Xe
%     dX(4) = -D_Xe / (rho*V); %Xe; %mantle Xe concentration
%     dX(5) =  D_Xe / (rho*V);
%     
%     if X(4)<0 && dX(4)<0
%         dX(4)=0;
%         X(4)=0;
%     end
%     
%     % Atmospheric Xe
%     dX(6) =  R_XesAtm / (rho*V); % atmospheric Xe in mantle
%     dX(7) = -R_XesAtm / (rho*V); % atmospheric Xe at surface
%     
%     if X(7)<0 && dX(7)<0
%         dX(7)=0;
%         X(7)=0;
%     end

% Xe      = X(4); % primordial Xe in mantle;  
% Xes     = X(5); % primordial Xe at surface -> ATM 128Xe/130Xe
% XeAtm   = X(6); % atmospheric Xe in mantle -> 128Xe in mantle
% XesAtm  = X(7); % atmospheric Xe at surface
  
%     % 130Xe DXecc
%     dX(4) = (R_XesAtm-D_Xe) / (rho*V) - DXecc*UrCC*exp(lamUr*(t_pd-t)/yr_s)*X(4)/Urt...
%         *(1/1e9/yr_s*( (t<time_cc)*0+(t>=time_cc)... %net CC growth
%         .*( ( 2*kcc*exp(kcc*(t-time_cc)/yr_s/1e9)...
%         .*(exp(kcc*(t_pd-time_cc)/yr_s/1e9)+1)...
%         ./((exp(kcc*(t_pd-time_cc)/yr_s/1e9)-1)*((exp(-kcc*(time_cc-t)/yr_s/1e9)+1)).^2) ) ...
%         +( Rccs/mcc+ (Rccp-Rccs)/mcc.*(1-exp(-krcc*(t-time_cc)/yr_s/1e9))./(1-exp(-krcc*(t_pd-time_cc)/yr_s/1e9)) )  ) ) )... %CC recycling
%         *mcc/V/rho; % actually dx/dt
%     dX(7) = 0;
    
%      dX(4) = (R_XesAtm-D_Xe) / (rho*V) - DXecc*X(4)*max( (UrCC*exp(lamUr*(t_pd-t)/yr_s)/Urt...
%         *(1/1e9/yr_s*( (t<time_cc)*0+(t>=time_cc)... %net CC growth
%         .*( ( 2*kcc*exp(kcc*(t-time_cc)/yr_s/1e9)...
%         .*(exp(kcc*(t_pd-time_cc)/yr_s/1e9)+1)...
%         ./((exp(kcc*(t_pd-time_cc)/yr_s/1e9)-1)*((exp(-kcc*(time_cc-t)/yr_s/1e9)+1)).^2) ) ...
%         +( Rccs/mcc+ (Rccp-Rccs)/mcc.*(1-exp(-krcc*(t-time_cc)/yr_s/1e9))./(1-exp(-krcc*(t_pd-time_cc)/yr_s/1e9)) )  ) ) )... %CC recycling
%         *mcc) - D0_mor,0 ) /V/rho; % actually dx/dt
%      dX(7) = 0;


    % 130Xe Dxecc=max(Dxecc0-Dxe,0)
    D_Xe = D0_mor*X(4);
    dX(4) = (R_XesAtm-D_Xe) / (rho*V) - DXecc*X(4)*D0_cc/V/rho - DXep*X(4)*D0_p/V/rho; % actually dx/dt
    dX(7) = 0;

    % 128Xe
    % ATM 128Xe/130Xe
    if t<2.6*1e9*yr_s
        dX(5)=-Xerpa*(130-128)*39/1000*1/(2.6*1e9*yr_s);
    else
        dX(5)=0;
    end    
    % 128Xe in mantle
    dX(6)=(R_XesAtm*Xes-D_Xe*XeAtm/Xe) / (rho*V) - DXecc*X(6)*D0_cc/V/rho - DXep*X(6)*D0_p/V/rho; 

    %Pu&Ur

    % cc model1: linear, extraction continues to today
%     dX(16)=-lamPu/yr_s*X(16)-UrCC*exp(lamUr*(t_pd-t)/yr_s)*Put/Urt*((t<(300e6*yr_s))*0+(t>=(300e6*yr_s))*(1/4.3e9/yr_s))*mcc/V/rho;
%     dX(17)=-lamUr/yr_s*X(17)-UrCC*exp(lamUr*(t_pd-t)/yr_s)*((t<(300e6*yr_s))*0+(t>=(300e6*yr_s))*(1/4.3e9/yr_s))*mcc/V/rho;

    % cc model2:sigmoid, extraction finished at 1Ga/2Ga/3Ga before
    dX(16)=-lamPu/yr_s*X(16)-UrCC*exp(lamUr*(t_pd-t)/yr_s)*Put/Urt...
        *NG_cc/V/rho;

    dX(17)=-lamUr/yr_s*X(17)-UrCC*exp(lamUr*(t_pd-t)/yr_s)...
        *NG_cc/V/rho;
    %



    % 136Xe
    % ATM 136Xe/130Xe
    if t<2.6*1e9*yr_s
        dX(8)=-Xe136rpa*(130-136)*39/1000*1/(2.6*1e9*yr_s);
    else
        dX(8)=0;
    end    
    % 136Xe in mantle
%     dX(9)=(R_XesAtm*Xe136r130atm-D_Xe*(Xe136+(lamPu/yr_s*YPu*X(16)+lamUr/yr_s*YUr*X(17))*par_timestep)/Xe) / (rho*V)+lamPu/yr_s*YPu*X(16)+lamUr/yr_s*YUr*X(17) ...
%         - DXecc*max( (UrCC*exp(lamUr*(t_pd-t)/yr_s)*X(9)/Urt...
%         *(1/1e9/yr_s*( (t<time_cc)*0+(t>=time_cc)... %net CC growth
%         .*( ( 2*kcc*exp(kcc*(t-time_cc)/yr_s/1e9)...
%         .*(exp(kcc*(t_pd-time_cc)/yr_s/1e9)+1)...
%         ./((exp(kcc*(t_pd-time_cc)/yr_s/1e9)-1)*((exp(-kcc*(time_cc-t)/yr_s/1e9)+1)).^2) ) ...
%         +( Rccs/mcc+ (Rccp-Rccs)/mcc.*(1-exp(-krcc*(t-time_cc)/yr_s/1e9))./(1-exp(-krcc*(t_pd-time_cc)/yr_s/1e9)) )  ) ) )... %CC recycling
%         *mcc)-D_Xe*(Xe136+(lamPu/yr_s*YPu*X(16)+lamUr/yr_s*YUr*X(17))*par_timestep)/Xe, 0)/V/rho; 
      
      dX(9)=( R_XesAtm*Xe136r130atm - ...
           ( D_Xe/Xe + DXecc*D0_cc + DXep*D0_p ) ...
          *( Xe136+(lamPu/yr_s*YPu*X(16)+lamUr/yr_s*YUr*X(17))*par_timestep) ) ...
          / (rho*V)...
          +lamPu/yr_s*YPu*X(16)+lamUr/yr_s*YUr*X(17); 



    % 131Xe
    % ATM 131Xe/130Xe
    if t<2.6*1e9*yr_s
        dX(10)=-Xe131rpa*(130-131)*39/1000*1/(2.6*1e9*yr_s);
    else
        dX(10)=0;
    end    
    % 131Xe in mantle
    dX(11)=( R_XesAtm*Xe131r130atm- ...
        ( D_Xe/Xe + DXecc*D0_cc + DXep*D0_p )...
        *(Xe131+(lamPu/yr_s*YPu131*X(16)+lamUr/yr_s*YUr131*X(17))*par_timestep)) ...
        / (rho*V) ...
        +lamPu/yr_s*YPu131*X(16)+lamUr/yr_s*YUr131*X(17); 

    % 132Xe
    % ATM 132Xe/130Xe
    if t<2.6*1e9*yr_s
        dX(12)=-Xe132rpa*(130-132)*39/1000*1/(2.6*1e9*yr_s);
    else
        dX(12)=0;
    end    
    % 132Xe in mantle
    dX(13)=( R_XesAtm*Xe132r130atm - ...
        ( D_Xe/Xe + DXecc*D0_cc + DXep*D0_p)...
        *(Xe132+(lamPu/yr_s*YPu132*X(16)+lamUr/yr_s*YUr132*X(17))*par_timestep) ) ...
        / (rho*V) ...
        +lamPu/yr_s*YPu132*X(16)+lamUr/yr_s*YUr132*X(17); 

    % 134Xe
    % ATM 134Xe/130Xe
    if t<2.6*1e9*yr_s
        dX(14)=-Xe134rpa*(130-134)*39/1000*1/(2.6*1e9*yr_s);
    else
        dX(14)=0;
    end    
    % 134Xe in mantle
    dX(15)=( R_XesAtm*Xe134r130atm-...
         ( D_Xe/Xe + DXecc*D0_cc + DXep*D0_p)...
        *( Xe134+(lamPu/yr_s*YPu134*X(16)+lamUr/yr_s*YUr134*X(17))*par_timestep) ) ...
        / (rho*V)...
        +lamPu/yr_s*YPu134*X(16)+lamUr/yr_s*YUr134*X(17); 

    
else
end


% if Fd_sw
%     % if X(3)<=1e-3 && dX(3)<0
%     if X(3)<=(1) && dX(3)<0
%         % if Ms < 1e-3 AND the rate of change of the surface ocean mass is -ve
%         X(3)  = 0; % if Ms is -ve, make it zero i.e. no water in the surface
%         dX(3) = (1/(M_oc*1e6))*(D); % surface ocean can increase due to degassing
%         dX(1) = -D./(rho.*V); % mantle can decreases due to degassing
%         % but as no water at the surface, no water can be subducted into the
%         % mantle
%         if (length(X)>3)
%             dX(4)=-D_Xe / (rho*V);%130Xe
%             dX(6)=-D_Xe*XeAtm/Xe / (rho*V);%128Xe
%             dX(9)=-D_Xe*(Xe136+lamPu/yr_s*YPu*X(16)+lamUr/yr_s*YUr*X(17))/Xe / (rho*V)+lamPu/yr_s*YPu*X(16)+lamUr/yr_s*YUr*X(17);%136Xe
%             dX(11)=-D_Xe*(Xe131+lamPu/yr_s*YPu131*X(16)+lamUr/yr_s*YUr131*X(17))/Xe / (rho*V)+lamPu/yr_s*YPu131*X(16)+lamUr/yr_s*YUr131*X(17);%131Xe
%             dX(13)=-D_Xe*(Xe132+lamPu/yr_s*YPu132*X(16)+lamUr/yr_s*YUr132*X(17))/Xe / (rho*V)+lamPu/yr_s*YPu132*X(16)+lamUr/yr_s*YUr132*X(17);%132Xe
%             dX(15)=-D_Xe*(Xe134+lamPu/yr_s*YPu134*X(16)+lamUr/yr_s*YUr134*X(17))/Xe / (rho*V)+lamPu/yr_s*YPu134*X(16)+lamUr/yr_s*YUr134*X(17);%134Xe
%         end
%     end
% 
% 
%     if X(1)<=0.1 && dX(1)<0
%     % if X(1)<=(1e6 * (Ms_init-1).*1.39e21/(3500*V)) && dX(1)<0
%         % if Xm < 0.1 AND the rate of change of the smantle water content is -ve
%         X(1)  = 0; % if X is -ve make it zero i.e. no water in the mantle 
%         dX(3) = (1/(M_oc*1e6))*(-R); % surface ocean can decrease due to regassing
%         dX(1) = R./(rho.*V); % mantle can increase due to regassing
%         % make both the change in X and Ms = 0
%     end
% 
% 
% else
    if X(3)<=1e-3 && dX(3)<0
        % if Ms < 1e-3 AND the rate of change of the surface ocean mass is -ve
        dX(3)=0;
        dX(1)=0;
         X(3)=0;
        % make both the change in X and Ms = 0
        % cause at this time, R is not zero, but also mantle rocks (R=D)
        % for Xe, R is propotional to R_water, in which Xp=X(1) now
        if (length(X)>3)
            dX(4)=(R_XesAtm/Xp*min(X(1),Xp)-D_Xe) / (rho*V) - DXecc*X(4)*D0_cc/V/rho - DXep*X(4)*D0_p/V/rho;

            dX(6)=(R_XesAtm/Xp*min(X(1),Xp)*Xes-D_Xe*XeAtm/Xe) / (rho*V)- DXecc*X(6)*D0_cc/V/rho - DXep*X(6)*D0_p/V/rho; 

            %dX(9)=0;%136Xe
            dX(9)=(R_XesAtm/Xp*min(X(1),Xp)*Xe136r130atm-D_Xe*(Xe136+lamPu/yr_s*YPu*X(16)+lamUr/yr_s*YUr*X(17))/Xe) / (rho*V)+lamPu/yr_s*YPu*X(16)+lamUr/yr_s*YUr*X(17) ...
                - DXecc*X(9)*D0_cc/V/rho - DXep*X(9)*D0_p/V/rho; 

            %dX(11)=0;%131Xe
            dX(11)=(R_XesAtm/Xp*min(X(1),Xp)*Xe131r130atm-D_Xe*(Xe131+lamPu/yr_s*YPu131*X(16)+lamUr/yr_s*YUr131*X(17))/Xe) / (rho*V)+lamPu/yr_s*YPu131*X(16)+lamUr/yr_s*YUr131*X(17) ...
                  - DXecc*X(11)*D0_cc/V/rho - DXep*X(11)*D0_p/V/rho; 
            
            %dX(13)=0;%132Xe
            dX(13)=(R_XesAtm/Xp*min(X(1),Xp)*Xe132r130atm-D_Xe*(Xe132+lamPu/yr_s*YPu132*X(16)+lamUr/yr_s*YUr132*X(17))/Xe) / (rho*V)+lamPu/yr_s*YPu132*X(16)+lamUr/yr_s*YUr132*X(17) ...
                - DXecc*X(13)*D0_cc/V/rho - DXep*X(13)*D0_p/V/rho; 
            

            %dX(15)=0;%134Xe
            dX(15)=(R_XesAtm/Xp*min(X(1),Xp)*Xe134r130atm-D_Xe*(Xe134+lamPu/yr_s*YPu134*X(16)+lamUr/yr_s*YUr134*X(17))/Xe) / (rho*V)+lamPu/yr_s*YPu134*X(16)+lamUr/yr_s*YUr134*X(17) ...
                  - DXecc*X(15)*D0_cc/V/rho - DXep*X(15)*D0_p/V/rho; 
        end
    end
    % 
    if X(1)<=0.1 && dX(1)<0
        % if Xm < 0.1 AND the rate of change of the smantle water content is -ve
        dX(3)=0;
        dX(1)=0;
        X(1)=0;
        % make both the change in X and Ms = 0
        % ? but at this time D!=R?
        % disp('a');
    end
    
    if (length(X)>3)
       if X(4)<=0.1 && dX(4)<0
            dX(4)=0;
            X(4)=0;
       end
       if X(6)<=0.1 && dX(6)<0
            dX(6)=0;
            X(6)=0;
       end
       if X(9)<=0.1 && dX(9)<0
            dX(9)=0;
            X(9)=0;
       end
       if X(11)<=0.1 && dX(11)<0
            dX(11)=0;
            X(11)=0;
       end
       if X(13)<=0.1 && dX(13)<0
            dX(13)=0;
            X(13)=0;
       end
       if X(15)<=0.1 && dX(15)<0
            dX(15)=0;
            X(15)=0;
       end
    end

%end





dX=dX';
%{
% if Fd_sw
%     % This will switch off degassing if there is less than 4 oceans in the
%     % mantle
%     if X(1)<1755 && dX(1)<0
%         X(1) = 1755; 
%         dX(1) = R/(rho*V); %Xmdot = kg/s
%         dX(3) = (1/(M_oc*1e6))*(-R);
%     end
% 
%     % This will switch off regassing if there is less than 4 oceans in the
%      mantle
%     if X(3)<1 && dX(3)<0
%         X(3) = 1; 
%         dX(3) = (1/(M_oc*1e6))*(D);
%         dX(1) = -D/(rho*V); %Xmdot = kg/s
%     end
% 
% end

% % This will swtich off degassing if there is no water in the mantle
% if X(1)<0.1 && dX(1)<0
%     X(1) = 0; 
%     dX(1) = R/(rho*V); %Xmdot = kg/s
%     dX(3) = (1/(M_oc*1e6))*(-R);
% end
% 
% % This will swtich off regassing if there is no water at the suface 1e-3
% if X(3)<1e-3 && dX(3)<0
%     X(3) = 0; 
%     dX(3) = (1/(M_oc*1e6))*(D);
%     dX(1) = -D/(rho*V); %Xmdot = kg/s
% end

%     A = 1e-10;
%     dX(4) = C * D * 1e-3 ; % = C * (D)/(rho*V) * (rho*V) * 1e-6 * 1e3; % g
%     dX(4) = - A * D/(rho*V) * 18e-6; %Xe; %mantle Xe concentration
%}