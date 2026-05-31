
par_solver_new;
tmax = t_pd/yr_s/1e9;% the age of solar system, in unit Ga

%% % initial CC
dt = 0.001;% length of each timestep
t = 0:dt:tmax;
nt = length(t);% number of timesteps
t = t';

load('cc/cc_output_for_figure.mat');

NGcc_cc = cc_output.NG_cc;
Gcc_cc = cc_output.G_cc;
Rcc_cc = cc_output.R_cc;
F_cc = cc_output.F_cc;
S_cc = cc_output.S_cc;


%%
[nt,size_F]=size(NGcc_cc);
tic;
[cc_output_margin.NGcc_cc_5,cc_output_margin.NGcc_cc_25,...
    cc_output_margin.NGcc_cc_50,cc_output_margin.NGcc_cc_75,cc_output_margin.NGcc_cc_95] = calculate_percentile_fun(NGcc_cc,nt,t,size_F);
[cc_output_margin.Gcc_cc_5,cc_output_margin.Gcc_cc_25,...
    cc_output_margin.Gcc_cc_50,cc_output_margin.Gcc_cc_75,cc_output_margin.Gcc_cc_95] = calculate_percentile_fun(Gcc_cc,nt,t,size_F);
[cc_output_margin.Rcc_cc_5,cc_output_margin.Rcc_cc_25,...
    cc_output_margin.Rcc_cc_50,cc_output_margin.Rcc_cc_75,cc_output_margin.Rcc_cc_95] = calculate_percentile_fun(Rcc_cc,nt,t,size_F);
[cc_output_margin.F_cc_5,cc_output_margin.F_cc_25,...
    cc_output_margin.F_cc_50,cc_output_margin.F_cc_75,cc_output_margin.F_cc_95] = calculate_percentile_fun(F_cc,nt,t,size_F);
[cc_output_margin.S_cc_5,cc_output_margin.S_cc_25,...
    cc_output_margin.S_cc_50,cc_output_margin.S_cc_75,cc_output_margin.S_cc_95] = calculate_percentile_fun(S_cc,nt,t,size_F);
elapsedTime = toc;  % s
disp(['used time = ', num2str(elapsedTime), ' s']);

[cc_output_margin.NGcc_lower,cc_output_margin.NGcc_upper] = calculate_margin_fun(NGcc_cc,nt,t,size_F);
[cc_output_margin.Gcc_lower,cc_output_margin.Gcc_upper] = calculate_margin_fun(Gcc_cc,nt,t,size_F);
[cc_output_margin.Rcc_lower,cc_output_margin.Rcc_upper] = calculate_margin_fun(Rcc_cc,nt,t,size_F);
[cc_output_margin.F_cc_lower,cc_output_margin.F_cc_upper] = calculate_margin_fun(F_cc,nt,t,size_F);
[cc_output_margin.S_cc_lower,cc_output_margin.S_cc_upper] = calculate_margin_fun(S_cc,nt,t,size_F);
elapsedTime = toc;  % s
disp(['used time = ', num2str(elapsedTime), ' s']);



tic;
[nt,size_Xe]=size(F_Xe);
[cc_output_margin.NGcc0_Xe_5,cc_output_margin.NGcc0_Xe_25,...
    cc_output_margin.NGcc0_Xe_50,cc_output_margin.NGcc0_Xe_75,cc_output_margin.NGcc0_Xe_95] = calculate_percentile_fun(NGcc0_Xe,nt,t,size_Xe);
[cc_output_margin.Gcc0_Xe_5,cc_output_margin.Gcc0_Xe_25,...
    cc_output_margin.Gcc0_Xe_50,cc_output_margin.Gcc0_Xe_75,cc_output_margin.Gcc0_Xe_95] = calculate_percentile_fun(Gcc0_Xe,nt,t,size_Xe);
[cc_output_margin.Rcc0_Xe_5,cc_output_margin.Rcc0_Xe_25,...
    cc_output_margin.Rcc0_Xe_50,cc_output_margin.Rcc0_Xe_75,cc_output_margin.Rcc0_Xe_95] = calculate_percentile_fun(Rcc0_Xe,nt,t,size_Xe);
[cc_output_margin.F_Xe_5,cc_output_margin.F_Xe_25,...
    cc_output_margin.F_Xe_50,cc_output_margin.F_Xe_75,cc_output_margin.F_Xe_95] = calculate_percentile_fun(F_Xe,nt,t,size_Xe);
[cc_output_margin.S_Xe_5,cc_output_margin.S_Xe_25,...
    cc_output_margin.S_Xe_50,cc_output_margin.S_Xe_75,cc_output_margin.S_Xe_95] = calculate_percentile_fun(S_Xe,nt,t,size_Xe);
elapsedTime = toc;  % s
disp(['used time = ', num2str(elapsedTime), ' s']);


[cc_output_margin.NGcc_Xe_lower,cc_output_margin.NGcc_Xe_upper] = calculate_margin_fun(NGcc0_Xe,nt,t,size_Xe);
[cc_output_margin.Gcc_Xe_lower,cc_output_margin.Gcc_Xe_upper] = calculate_margin_fun(Gcc0_Xe,nt,t,size_Xe);
[cc_output_margin.Rcc_Xe_lower,cc_output_margin.Rcc_Xe_upper] = calculate_margin_fun(Rcc0_Xe,nt,t,size_Xe);
[cc_output_margin.F_Xe_lower,cc_output_margin.F_Xe_upper] = calculate_margin_fun(F_Xe,nt,t,size_Xe);
[cc_output_margin.S_Xe_lower,cc_output_margin.S_Xe_upper] = calculate_margin_fun(S_Xe,nt,t,size_Xe);
elapsedTime = toc;  % s
disp(['used time = ', num2str(elapsedTime), ' s']);


% save('cc/cc_output_for_figure_margin.mat','cc_output_margin','-v7.3');


%% % initial T
dt = 0.001;% length of each timestep
t = 0:dt:tmax;
nt = length(t);% number of timesteps
t = t';

load('thermal/thermal_output_for_figure.mat');

% load
NGcc_T=T_output.NG_cc  ;
Gcc_T=T_output.G_cc  ;
Rcc_T=T_output.R_cc  ;
% F_T=T_output.F_T  ;
% S_T=T_output.S_T  ;
Tp_T=T_output.Tp_T;
Tm_T=T_output.Tm_T;
Qs_T=T_output.Qs_T;
Qc_T=T_output.Qc_T;
H_T=T_output.H_T;
Ra_T=T_output.Ra_T;

[nt,size_T]=size(Tp_T);

% Calculate the boundary of evolution curves
tic;
% [NGcc_T_5,NGcc_T_25,NGcc_T_50,NGcc_T_75,NGcc_T_95] = calculate_percentile_fun(NGcc_T,nt,t,size_T);
% [Rcc_T_5,Rcc_T_25,Rcc_T_50,Rcc_T_75,Rcc_T_95] = calculate_percentile_fun(Rcc_T,nt,t,size_T);
% [Gcc_T_5,Gcc_T_25,Gcc_T_50,Gcc_T_75,Gcc_T_95] = calculate_percentile_fun(Gcc_T,nt,t,size_T);
[thermal_output_margin.Tp_T_5,thermal_output_margin.Tp_T_25,...
    thermal_output_margin.Tp_T_50,thermal_output_margin.Tp_T_75,thermal_output_margin.Tp_T_95] = calculate_percentile_fun(Tp_T,nt,t,size_T);
[thermal_output_margin.Qs_T_5,thermal_output_margin.Qs_T_25,...
    thermal_output_margin.Qs_T_50,thermal_output_margin.Qs_T_75,thermal_output_margin.Qs_T_95] = calculate_percentile_fun(Qs_T/1e12,nt,t,size_T);
elapsedTime = toc;  % s
disp(['used time = ', num2str(elapsedTime), ' s']);


[thermal_output_margin.Tp_T_lower,thermal_output_margin.Tp_T_upper]=calculate_margin_fun(Tp_T,nt,t,size_T);
[thermal_output_margin.Qs_T_lower,thermal_output_margin.Qs_T_upper]=calculate_margin_fun(Qs_T,nt,t,size_T);
[thermal_output_margin.Qc_T_lower,thermal_output_margin.Qc_T_upper]=calculate_margin_fun(Qc_T,nt,t,size_T);
[thermal_output_margin.H_T_lower,thermal_output_margin.H_T_upper]=calculate_margin_fun(H_T,nt,t,size_T);
[thermal_output_margin.Ra_T_lower,thermal_output_margin.Ra_T_upper]=calculate_margin_fun(Ra_T,nt,t,size_T);
elapsedTime = toc;  % s
disp(['used time = ', num2str(elapsedTime), ' s']);

% save('thermal/thermal_output_for_figure_margin.mat','thermal_output_margin','-v7.3');


%% % after Xe
default_input_for_degassing;
t=Simulation_input.time_series;
t=t/yr_s/1e9;% to Ga
nt = length(t);
ntd=nt;

% load('Xe/Xe_output_for_figure.mat');

% save
% cc_Xe=Xe_output.cc_Xe  ;
% NGcc_Xe=Xe_output.NG_cc  ;
% Gcc_Xe=Xe_output.G_cc  ;
% Rcc_Xe=Xe_output.R_cc  ;
% % F_Xe=Xe_output.F_Xe  ;
% % S_Xe=Xe_output.S_Xe  ;
% Tm_Xe=Xe_output.Tm_Xe;
% Qs_Xe=Xe_output.Qs_Xe;
% Qc_Xe=Xe_output.Qc_Xe;
% H_Xe=Xe_output.H_Xe;
% Ra_Xe=Xe_output.Ra_Xe;
% U_Xe=Xe_output.U_Xe ; % velocity
% 
% Xe=Xe_output.Xe;
% Xe128=Xe_output.Xe128;
% Xe131=Xe_output.Xe131;
% Xe132=Xe_output.Xe132;
% Xe134=Xe_output.Xe134;
% Xe136=Xe_output.Xe136;
% Pu=Xe_output.Pu;
% Ur=Xe_output.Ur;

[nt,size_Xe]=size(Tp_Xe);
ntd=nt;

%% calculate the boundary of evoltion curves
tic;
[Xe_output_margin.cc_Xe_5,Xe_output_margin.cc_Xe_25,...
    Xe_output_margin.cc_Xe_50,Xe_output_margin.cc_Xe_75,Xe_output_margin.cc_Xe_95] = calculate_percentile_fun(cc_Xe,nt,t,size_Xe);
[Xe_output_margin.NGcc_Xe_5,Xe_output_margin.NGcc_Xe_25,...
    Xe_output_margin.NGcc_Xe_50,Xe_output_margin.NGcc_Xe_75,Xe_output_margin.NGcc_Xe_95] = calculate_percentile_fun(NGcc_Xe,nt,t,size_Xe);
[Xe_output_margin.Rcc_Xe_5,Xe_output_margin.Rcc_Xe_25,...
    Xe_output_margin.Rcc_Xe_50,Xe_output_margin.Rcc_Xe_75,Xe_output_margin.Rcc_Xe_95] = calculate_percentile_fun(Rcc_Xe,nt,t,size_Xe);
% [Xe_output_margin.F_Xe_5,Xe_output_margin.F_Xe_25,...
%     Xe_output_margin.F_Xe_50,Xe_output_margin.F_Xe_75,Xe_output_margin.F_Xe_95] = calculate_percentile_fun(F_Xe,nt,t,size_Xe);
% [Xe_output_margin.S_Xe_5,Xe_output_margin.S_Xe_25,...
%     Xe_output_margin.S_Xe_50,Xe_output_margin.S_Xe_75,Xe_output_margin.S_Xe_95] = calculate_percentile_fun(S_Xe,nt,t,size_Xe);

[Xe_output_margin.Tp_Xe_5,Xe_output_margin.Tp_Xe_25,...
    Xe_output_margin.Tp_Xe_50,Xe_output_margin.Tp_Xe_75,Xe_output_margin.Tp_Xe_95] = calculate_percentile_fun(Tp_Xe,ntd,t,size_Xe);
% [Qs_Xe_5,Qs_Xe_25,Qs_Xe_50,Qs_Xe_75,Qs_Xe_95] = calculate_percentile_fun(Qs_Xe/1e12,ntd,t,size_Xe);

[Xe_output_margin.D0_mor_5,Xe_output_margin.D0_mor_25,...
    Xe_output_margin.D0_mor_50,Xe_output_margin.D0_mor_75,Xe_output_margin.D0_mor_95] = calculate_percentile_fun(D0_mor,ntd,t,size_Xe);
[Xe_output_margin.D0_cc_5,Xe_output_margin.D0_cc_25,...
    Xe_output_margin.D0_cc_50,Xe_output_margin.D0_cc_75,Xe_output_margin.D0_cc_95] = calculate_percentile_fun(D0_cc,ntd,t,size_Xe);
[Xe_output_margin.D0_p_5,Xe_output_margin.D0_p_25,...
    Xe_output_margin.D0_p_50,Xe_output_margin.D0_p_75,Xe_output_margin.D0_p_95] = calculate_percentile_fun(D0_p,ntd,t,size_Xe);


[Xe_output_margin.Xe_5,Xe_output_margin.Xe_25,...
    Xe_output_margin.Xe_50,Xe_output_margin.Xe_75,Xe_output_margin.Xe_95] = calculate_percentile_fun(Xe,ntd,t,size_Xe);
[Xe_output_margin.Xe128_5,Xe_output_margin.Xe128_25,...
    Xe_output_margin.Xe128_50,Xe_output_margin.Xe128_75,Xe_output_margin.Xe128_95] = calculate_percentile_fun(Xe128,ntd,t,size_Xe);
[Xe_output_margin.Xe136_5,Xe_output_margin.Xe136_25,...
    Xe_output_margin.Xe136_50,Xe_output_margin.Xe136_75,Xe_output_margin.Xe136_95] = calculate_percentile_fun(Xe136,ntd,t,size_Xe);

[Xe_output_margin.Tp_Xe_lower,Xe_output_margin.Tp_Xe_upper] = calculate_margin_fun(Tp_Xe,ntd,t,size_Xe);
[Xe_output_margin.Xe_lower,Xe_output_margin.Xe_upper] = calculate_margin_fun(Xe,ntd,t,size_Xe);

[Xe_output_margin.D0_mor_lower,Xe_output_margin.D0_mor_upper] = calculate_margin_fun(D0_mor,ntd,t,size_Xe);
[Xe_output_margin.D0_cc_lower,Xe_output_margin.D0_cc_upper] = calculate_margin_fun(D0_cc,ntd,t,size_Xe);
[Xe_output_margin.D0_p_lower,Xe_output_margin.D0_p_upper] = calculate_margin_fun(D0_p,ntd,t,size_Xe);

[Xe_output_margin.Pu_lower,Xe_output_margin.Pu_upper] = calculate_margin_fun(Pu,ntd,t,size_Xe);
[Xe_output_margin.Ur_lower,Xe_output_margin.Ur_upper] = calculate_margin_fun(Ur,ntd,t,size_Xe);
elapsedTime = toc;  % s
disp(['used time = ', num2str(elapsedTime), ' s']);

% save('Xe/Xe_output_for_figure_margin1.mat','Xe_output_margin','-v7.3');




