%% Parameters

yr_s=365.25*24*60*60; %s in a year
t_pd = 4.6e9*yr_s;
par_steps=4000;
par_timestep=0;%t_pd/par_steps;
Ts=273;%surface temperature [K]

% core heat flux
% Qc_pd = 15e12; % unit: W
% dQc = 2e12; % unit: W

% D - Degassing
% A factor - step change of degassing  or regassing
DXecc=1;
Dwcc=1;
Dwp=1;
DXep=1;
Rxe=0.9;
ADt=4.0e9;
ADf=1;
ARt=0;%3e9
ARf1=1;%0
ARf2=1;%4.6e9/(4.6e9-ARt)

%maximum plate age for Labrosse scaling
age_max=180e6*yr_s;%180Myr


% Planet Size
r_p = 6371e3; %m radius of planet
d = 2890e3;     %m depth of mantle
r_c = r_p-d; %m^3 radius of planet core
V = (4*pi/3)*(r_p^3-r_c^3);   %m^3 mantle volume
M = 4e24; %kg

% Viscosity - Mei and Kohlstedt (2000a,b)
r_eta = 1; n_eta = 1; San_r_eta=1.2;
eta0 = 1e21; %Pas viscosity
%eta00=2e15;
San_eta0 = 5.2e18; %1.7e17; %dimensionless calibration constant
E = 300e3;      %J/mol activation energy 
McG_E = 581.98e3; %McGovern and Schubert 1980
San_E = 480e3;
Rg = 8.314;      %J/mol.K ideal gas constant


% Rayleigh Number
alpha_ra = 2.5e-5;% 2.5e-5; %1/K thermal expansion coefficient
g = 9.8;       %m/s^2 gravitational acceleration
%rho = 3500;     %kg/m^3 density
rho = 4400; %kg/m^3 density used to calculate mantle mass (as a bulk to be mixed by D&R)
rhom = 3500; %kg/m^3 (upper mantle) density used to calculate degassing & regassing mass flux 
rhoa = 3500;  %kg/m^3 density used to scale Ra, k, eta/rhoa

% plume degassing flux
fpb_pd = 50e3; %plume bouyancy kg/s ref King (2014)
pdT = 300; % plume thermal anomaly K Steinberger (2024, "Why Are Plume Excess Temperatures Much Less …”)
D0p_pd = fpb_pd/(alpha_ra*pdT);

% Plate Velocity
%L = 2*d; %m average plate length/ average convection cell length
% actually not affect the simulation cause tao_max=L/velocity has no L; 
% Take pi to match present-day tao_max~180-200 Myr, velocity~5cm/yr 
L = pi*d; 


% Depth of Melting - Hirschmann et al. (2009)
% constants
z1 = 286;    %m/K
z2 = 164;    %m/ppm
z3 = -3.266e5;   %m




% Qs - Total Surface Heatflow
S = 0.71*4*pi*r_p^2; %m^2 planet surface area (oceanic part)
Scc=0.29*4*pi*r_p^2;
Kc = 3; 
Cp = 1000;
k = Kc/(Cp*rhoa);     %m2/s thermal diffusivity
Ra_cr=1100;




% Water Fugacity
% experimentally determined constants
c0 = -7.9859; c1 = 4.3559; c2 = -0.5742; c3 = 0.0337;
M_oc = 1.39e21; %kg mass of ocean today

% scaling constant for viscosity

% H_Si=100/6; %conversion from ppm to H/10^6Si atoms for Li et al.(2008) equation
% A_Coh = 500*H_Si; %500ppm to H/10^6Si
% A_F_h2o = exp(c0 + (c1*log(A_Coh)) + (c2*(log(A_Coh))^2) + (c3*(log(A_Coh))^3)); %MPa to Pa
% A_eta_w = (A_F_h2o^(-r_eta)); %1/Pa
% A_eta_T = exp(E/(Rg*n_eta*1573)); %1300C
% Acre = eta0/abs(A_eta_w*A_eta_T); %scaling constant for viscosity: Pas.Pa

H_Si = 100/6; %conversion from ppm to H/10^6Si atoms for Li et al.(2008) equation
A_Coh = 500*H_Si; %500ppm to H/10^6Si
A_F_h2o = exp(c0 + (c1*log(A_Coh)) + (c2*(log(A_Coh))^2) + (c3*(log(A_Coh))^3)); %MPa to Pa
A_eta_w = (A_F_h2o^(-r_eta)); %1/Pa
A_eta_T = exp(E/(Rg*n_eta*1573)); %1300C
A = eta0/abs(A_eta_w*A_eta_T); %scaling constant for viscosity: Pas.Pa

% uncorrected code
H_Si_un=6.3e-4; %conversion from ppm to H/10^6Si atoms for Li et al.(2008) equation
A_Coh_un = 500*H_Si_un; %500ppm to H/10^6Si
A_F_h2o_un = 1e6*exp(c0 + (c1*log(A_Coh_un)) + (c2*(log(A_Coh_un))^2) + (c3*(log(A_Coh_un))^3)); %MPa to Pa
A_eta_w_un = (A_F_h2o_un^(-r_eta)); %1/Pa
A_eta_T_un = exp(E/(Rg*n_eta*1300)); %1300C
A_un = eta0/abs(A_eta_w_un*A_eta_T_un); %scaling constant for viscosity: Pas.Pa

% uncorrected code corrected T
H_Si_unT=6.3e-4; %conversion from ppm to H/10^6Si atoms for Li et al.(2008) equation
A_Coh_unT = 500*H_Si_un; %500ppm to H/10^6Si
A_F_h2o_unT = 1e6*exp(c0 + (c1*log(A_Coh_unT)) + (c2*(log(A_Coh_unT))^2) + (c3*(log(A_Coh_unT))^3)); %MPa to Pa
A_eta_w_unT = (A_F_h2o_unT^(-r_eta)); %1/Pa
A_eta_T_unT = exp(E/(Rg*n_eta*1573)); %1300C
A_unT = eta0/abs(A_eta_w_unT*A_eta_T_unT); %scaling constant for viscosity: Pas.Pa

% sandu et al. (2011)
San_Acre=90; 

% water concentration in plate
serp_frac=0.2;
rho_serp=3100; %kg/m3
rho_dry=3300;  %kg/m3
Xdry=0;
Xserp=1.3e4;   %ppm ~ 1.3wt%
Xp = 2.056e+03; % Xp = (rho_serp*serp_frac*Xserp)./((rho_serp*serp_frac)+rho_dry) ppm
%Xp = 1.3e4;

NA=6.02e23;
% Xe isotope composition
Xercc=0.5073; % initial 128Xe/130Xe for mantle (ref cc = carbacenous chondrite)
Xe136rcc=1.988; % initial 136Xe/130Xe for mantle (ref cc)
Xe131rcc=5.043; % initial 131Xe/130Xe for mantle (ref cc)
Xe132rcc=6.150; % initial 132Xe/130Xe for mantle (ref cc)
Xe134rcc=2.359; % initial 134Xe/130Xe for mantle (ref cc)

Xeren=0.50812; % initial 128Xe/130Xe for mantle (ref en = enstatite Patzer & Schultz, 2002)
Xe136ren=1.92513; % initial 136Xe/130Xe for mantle (ref en)
Xe131ren=5.02835; % initial 131Xe/130Xe for mantle (ref en)
Xe132ren=6.1389; % initial 132Xe/130Xe for mantle (ref en)
Xe134ren=2.32357; % initial 134Xe/130Xe for mantle (ref en)


atmf=39/1000;
Xerpa=0.4715; % present day 128Xe/130Xe for ATM
Xe136rpa=2.176; % present day 136Xe/130Xe for ATM
Xe131rpa=5.213; % present day 131Xe/130Xe for ATM
Xe132rpa=6.607; % present day 132Xe/130Xe for ATM
Xe134rpa=2.563; % present day 134Xe/130Xe for ATM

% Pu and U para
%Urini=40e-9/238*NA;%21 ppb 238U for BSE unit [ppm] today-42.1958ppb initially;or 30.8-62.8756
UCC=1.3e-6; %1.3ppm mass ratio ref Rudnick and Gao, 2003
U238CC=UCC*0.993;
U235CC=UCC*0.007;
ThCC=5.6e-6;%5.6ppm
KCC=1.5e-2*1.19e-4;%1.5wt%K*1.19e-4 


%for fission, need use atoms/g unit
UrCC=U238CC/238*NA;%1.3ppm U in CC - change into atoms/g
mcc_val=2.2e22;%2.2e22;%mass (kg) of CC
YPu=7e-5;%136Xe
YUr=3.43e-8;%136Xe
YPu132=YPu/1.738;%132Xe
YUr132=YUr/1.120;%132Xe
YPu131=YPu132*0.1449;%131Xe
YUr131=YUr132*0.2777;%131Xe
YPu134=YPu132*1.437;%134Xe
YUr134=YUr132*1.041;%134Xe

apg_ppm=2.16e-16;%130/6.02e23/1e-6
%128/130: apg ratio * 128/130 = ppm ratio
%xe128_mreal=xe128*128/130;
%ratio_real=xe128_mreal/xe/128*130;










