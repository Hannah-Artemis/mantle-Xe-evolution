%% % calculate quatities based on successful params

%
ncores=4;
if isempty(gcp('nocreate'))
    parpool('local', ncores);  % 
end

%% set params
par_solver_new;
tmax = t_pd/yr_s/1e9;% the age of solar system, in unit Ga

Mcp = mcc_val;% mass of continental crust at present-day, in unit kg
para_Scale=2;% 2: Schubert
Msi=1.5;% 
SwXe = logical(1);
Xei=3.2e7; % atoms/g



%% % load data
success_F_struct=load('cc/out_cc_grwoth_success_lhs_MC_Krw.mat');
success_F=success_F_struct.success_F;
[size_F,nparam_cc]=size(success_F);

success_T_struct=load('thermal/out_thermal_success_lhs_test_fixK_U_MC_lowerT.mat');
success_T=success_T_struct.success_T;
[size_T,nparam_T]=size(success_T);

success_Xe_struct=load('Xe/out_Xe_success_lhs_new_3e7_test_fixK_U_MC_negativFD.mat');
success_Xe=success_Xe_struct.success_Xe;
[size_Xe,nparam_Xe]=size(success_Xe);

% find Krw for CC params
Krw_Xe=zeros(size_Xe,1);
for i=1:size_Xe
    idcc_mask=find(success_F(:,2)==success_Xe(i,2));
    Krw_Xe(i)=success_F(idcc_mask,1);
end

%% % initial CC
dt = 0.001;% length of each timestep
t = 0:dt:tmax;
nt = length(t);% number of timesteps
t = t';

%% Load in the observed formation & surface age distributions
% Formation age distribution data from Korenaga (2018a)
%data_formationage = load('korenaga18a_Tunmix_orig.dat');
data_formationage = readmatrix(...
    'cc/observation/cc_age.xlsx','Sheet','korenaga18a_Tunmix_orig');
% Surface age distribution data from Roberts & Spencer (2015)
data_zircon_surf = readmatrix(...
    'cc/observation/cc_age.xlsx','Sheet','korenaga18a_T_U_Pb');
% use load_FandS_fun function to set data in the same dimension as the time series
[F_Jun_same,S_Jun_same] = load_FandS_fun(t,nt,data_formationage,data_zircon_surf);

%%
% test running time: cc_grwoth_fun(0.12s) formation_surface_age(0.36s); 
% all on 6 cores : 1400s
% size_F=12;
tic;

futures1(size_F) = parallel.FevalFuture;
futures2(size_F) = parallel.FevalFuture;

F_cc = nan(nt,size_F);
S_cc = nan(nt,size_F);
NGcc_cc = nan(nt,size_F);
Rcc_cc = nan(nt,size_F);
Gcc_cc = nan(nt,size_F);
Krw = nan(nt,size_F);

% cc-curves
for l=1:size_F
    if mod(l,10)==0
        disp(['lcc10=' num2str(l) ' of ' num2str(size_F)]);% keep track of the calculation
    end
% independent main variables
    % Krw_factor_model = success_F(l,1);
    ts_model = success_F(l,1);
    kappa_g_model = success_F(l,2);
    kappa_r_model = success_F(l,3);    
    Rs_model = success_F(l,4);
    Rp_model = success_F(l,5);
    frw_model = success_F(l,6);
    % Calculate the dependent variable (crustal reworking rate)
    % Krw_s_model = Rs_model * Krw_factor_model;% initial Krw_factor
    % Calculate the corresponding crustal growth pattern
%     [NGcc_model,Rcc_model,Gcc_model,Krw_model] = CC_growth_fun2(t,ts_model,tmax,...
%     Mcp,kappa_g_model,Rp_model,Rs_model,kappa_r_model,frw_model);

    futures1(l) = parfeval(@CC_growth_fun2, 4, t,ts_model,tmax,...
     Mcp,kappa_g_model,Rp_model,Rs_model,kappa_r_model,frw_model);


    % Calculate formation age and surface age distributions and recalculate
    % crustal reworking rate using Formation_surface_age_fun function
% [F_model, S_model, m_tp,m,Krw_real_model,s] = Formation_surface_age_fun(t, Gcc_model, Rcc_model, NGcc_model, Krw_model);

%  futures2(l) = parfeval(@Formation_surface_age_fun, 2, t, Gcc_model,...
%      Rcc_model, NGcc_model, Krw_model);   


%     % save the results
%    F_cc(:,l) = F_model;
%    S_cc(:,l) = S_model;
%     NGcc_cc(:,l) = NGcc_model;
%     Rcc_cc(:,l) = Rcc_model;
%     Gcc_cc(:,l) = Gcc_model;
end % for l=1:itermax
elapsedTime = toc;  % s
disp(['used time = ', num2str(elapsedTime), ' s']);

% cc-curves
for i = 1:size_F
    if mod(i,10)==0
        disp(['lcc11=' num2str(i) ' of ' num2str(size_F)]);% keep track of the calculation
    end
    try
        %results(i, :) = fetchOutputs(futures(i));
        [NGcc_cc(:,i),Rcc_cc(:,i),Gcc_cc(:,i),Krw(:,i)]=fetchOutputs(futures1(i));
        % [F_cc(:,i),S_cc(:,i)]=fetchOutputs(futures2(i));
    catch ME
        warning('Error in task %d: %s', i, ME.message);
        % results(i, :) = nan(nt,1);  %
        NGcc_cc(:,i)=nan(nt,1);
        Rcc_cc(:,i)=nan(nt,1);
        Gcc_cc(:,i)=nan(nt,1);
        Krw(:,i)=nan(nt,1);
%         F_cc(:,i)=nan(nt,1);
%         S_cc(:,i)=nan(nt,1);
    end
end
elapsedTime = toc;  % s
disp(['used time = ', num2str(elapsedTime), ' s']);



%% % cc-age
for l=1:size_F
    if mod(l,10)==0
        disp(['lcc20=' num2str(l) ' of ' num2str(size_F)]);% keep track of the calculation
    end
% independent main variables
     NGcc_model = NGcc_cc(:,l);
     Rcc_model = Rcc_cc(:,l);
     Gcc_model = Gcc_cc(:,l);
     Krw_model = Krw(:,l);

    % Calculate formation age and surface age distributions and recalculate
    % crustal reworking rate using Formation_surface_age_fun function
 futures2(l) = parfeval(@Formation_surface_age_fun, 2, t, Gcc_model,...
     Rcc_model, NGcc_model, Krw_model);   
end % for l=1:itermax


for i = 1:size_F
    if mod(i,10)==0
        disp(['lcc21=' num2str(i) ' of ' num2str(size_F)]);% keep track of the calculation
    end
    try
        %results(i, :) = fetchOutputs(futures(i));
        % [NGcc_cc(:,i),Rcc_cc(:,i),Gcc_cc(:,i),Krw(:,i)]=fetchOutputs(futures1(i));
        [F_cc(:,i),S_cc(:,i)]=fetchOutputs(futures2(i));
    catch ME
        warning('Error in task %d: %s', i, ME.message);
        % results(i, :) = nan(nt,1);  %
%         NGcc_cc(:,i)=nan(nt,1);
%         Rcc_cc(:,i)=nan(nt,1);
%         Gcc_cc(:,i)=nan(nt,1);
%         Krw(:,i)=nan(nt,1);
         F_cc(:,i)=nan(nt,1);
         S_cc(:,i)=nan(nt,1);
    end
end
elapsedTime = toc;  % s
disp(['used time = ', num2str(elapsedTime), ' s']);

%% % Xe-cc
tic;
futures1(size_Xe) = parallel.FevalFuture;
futures2(size_Xe) = parallel.FevalFuture;

F_Xe = nan(nt,size_Xe);
S_Xe = nan(nt,size_Xe);
NGcc0_Xe = nan(nt,size_Xe);
Rcc0_Xe = nan(nt,size_Xe);
Gcc0_Xe = nan(nt,size_Xe);
Krw0_Xe = nan(nt,size_Xe);

% Xe-curves
for l=1:size_Xe
    if mod(l,10)==0
        disp(['lXecc10=' num2str(l) ' of ' num2str(size_Xe)]);% keep track of the calculation
    end
% independent main variables
    % Krw_factor_model = success_F(l,1);
    ts_model = success_Xe(l,1);
    kappa_g_model = success_Xe(l,2);
    kappa_r_model = success_Xe(l,3);    
    Rs_model = success_Xe(l,4);
    Rp_model = success_Xe(l,5);
    frw_model = Krw_Xe(l);
    % Calculate the dependent variable (crustal reworking rate)
    % Krw_s_model = Rs_model * Krw_factor_model;% initial Krw_factor
    % Calculate the corresponding crustal growth pattern
%     [NGcc_model,Rcc_model,Gcc_model,Krw_model] = CC_growth_fun2(t,ts_model,tmax,...
%     Mcp,kappa_g_model,Rp_model,Rs_model,kappa_r_model,frw_model);

    futures1(l) = parfeval(@CC_growth_fun2, 4, t,ts_model,tmax,...
     Mcp,kappa_g_model,Rp_model,Rs_model,kappa_r_model,frw_model);


end % for l=1:itermax
elapsedTime = toc;  % s
disp(['used time = ', num2str(elapsedTime), ' s']);


for i = 1:size_Xe
    if mod(i,10)==0
        disp(['lXecc11=' num2str(i) ' of ' num2str(size_Xe)]);% keep track of the calculation
    end
    try
        %results(i, :) = fetchOutputs(futures(i));
        [NGcc0_Xe(:,i),Rcc0_Xe(:,i),Gcc0_Xe(:,i),Krw0_Xe(:,i)]=fetchOutputs(futures1(i));
        % [F_cc(:,i),S_cc(:,i)]=fetchOutputs(futures2(i));
    catch ME
        warning('Error in task %d: %s', i, ME.message);
        % results(i, :) = nan(nt,1);  %
        NGcc0_Xe(:,i)=nan(nt,1);
        Rcc0_Xe(:,i)=nan(nt,1);
        Gcc0_Xe(:,i)=nan(nt,1);
        Krw0_Xe(:,i)=nan(nt,1);
%         F_cc(:,i)=nan(nt,1);
%         S_cc(:,i)=nan(nt,1);
    end
end
elapsedTime = toc;  % s
disp(['used time = ', num2str(elapsedTime), ' s']);



%% % age
for l=1:size_Xe
    if mod(l,10)==0
        disp(['lXecc20=' num2str(l) ' of ' num2str(size_Xe)]);% keep track of the calculation
    end
% independent main variables
     NGcc_model = NGcc0_Xe(:,l);
     Rcc_model = Rcc0_Xe(:,l);
     Gcc_model = Gcc0_Xe(:,l);
     Krw_model = Krw0_Xe(:,l);

    % Calculate formation age and surface age distributions and recalculate
    % crustal reworking rate using Formation_surface_age_fun function
 futures2(l) = parfeval(@Formation_surface_age_fun, 2, t, Gcc_model,...
     Rcc_model, NGcc_model, Krw_model);   
end % for l=1:itermax


for i = 1:size_Xe
    if mod(i,10)==0
        disp(['lXecc21=' num2str(i) ' of ' num2str(size_Xe)]);% keep track of the calculation
    end
    try
        %results(i, :) = fetchOutputs(futures(i));
        % [NGcc_cc(:,i),Rcc_cc(:,i),Gcc_cc(:,i),Krw(:,i)]=fetchOutputs(futures1(i));
        [F_Xe(:,i),S_Xe(:,i)]=fetchOutputs(futures2(i));
    catch ME
        warning('Error in task %d: %s', i, ME.message);
        % results(i, :) = nan(nt,1);  %
%         NGcc_cc(:,i)=nan(nt,1);
%         Rcc_cc(:,i)=nan(nt,1);
%         Gcc_cc(:,i)=nan(nt,1);
%         Krw(:,i)=nan(nt,1);
         F_Xe(:,i)=nan(nt,1);
         S_Xe(:,i)=nan(nt,1);
    end
end
elapsedTime = toc;  % s
disp(['used time = ', num2str(elapsedTime), ' s']);




%% % save
cc_output.NG_cc = NGcc_cc;
cc_output.G_cc = Gcc_cc;
cc_output.R_cc = Rcc_cc;
cc_output.F_cc = F_cc;
cc_output.S_cc = S_cc;

cc_output.NG_Xe = NGcc0_Xe;
cc_output.G_Xe = Gcc0_Xe;
cc_output.R_Xe = Rcc0_Xe;
cc_output.F_Xe = F_Xe;
cc_output.S_Xe = S_Xe;

% save('cc/cc_output_for_figure.mat','cc_output','-v7.3');
% save('cc/out_Xe_lhs_new_3e7_test_fixK_U_MC.mat', 'results', 'col_names','-v7.3');



%% % initial T
% all on 4 cores: 320s (6 cores will run out the storage)
eta_ref = 7.5e14;% fixed eta_ref from the relationship of present day Tp & Qs
para_Scale=2;% 2: Schubert
Hsf=2;
tmax = t_pd/yr_s/1e9;% the age of solar system, in unit Ga
dt = 0.001;% length of each timestep
t = 0:dt:tmax;
nt = length(t);% number of timesteps
t = t';

%% Load in the observed potential temperature history

% observation 1: present day potential temperature
Tp_pd=1350;
Tp_range=100;
Tp_obs1 = struct('Tp_pd', Tp_pd, 'Tp_range', Tp_range);

% observation 2: potential Temperature curve anchor points
data_Tp = readmatrix(...
    'thermal/observation/Tp.xlsx','Sheet','Herz data');
[Tp_anchorHerz1,Tp_anchorHerz2,Tp_anchorHerz3,Tp_anchorHerz4,...
    t_anchorHerz1,t_anchorHerz2,t_anchorHerz3,t_anchorHerz4,...
    t_Herz,Tp_Herz] = load_Tp_fun(data_Tp,tmax);

% observation 3: present day heat flux
 Qs_range=[26.7e12,34e12]; % W
% Qs_range=[33e12,45e12]; % W


% use anchor points only
Tp_obs2 = struct( ...
    'Tp_anchorHerz1', Tp_anchorHerz1, ...
    'Tp_anchorHerz2', Tp_anchorHerz2, ...
    'Tp_anchorHerz3', Tp_anchorHerz3, ...
    'Tp_anchorHerz4', Tp_anchorHerz4, ...
    't_anchorHerz1', t_anchorHerz1, ...
    't_anchorHerz2', t_anchorHerz2, ...
    't_anchorHerz3', t_anchorHerz3, ...
    't_anchorHerz4', t_anchorHerz4 ...
);


%%
tic;
% size_T=12;

futures1(size_T) = parallel.FevalFuture;
futures2(size_T) = parallel.FevalFuture;

NGcc_T = nan(nt,size_T);
Rcc_T = nan(nt,size_T);
Gcc_T = nan(nt,size_T);
Tm_T = nan(nt,size_T);
Tp_T = nan(nt,size_T);
Qs_T = nan(nt,size_T);
Qc_T = nan(nt,size_T);
H_T = nan(nt,size_T);
Ra_T = nan(nt,size_T);



for l=1:size_T
    if mod(l,10)==0
        disp(['lT10=' num2str(l) ' of ' num2str(size_T)]);% keep track of the calculation
    end
% independent main variables
    
    % p.time_cc_input,p.kappa_gcc,p.BSE,p.Qc_pd_input,p.dQc_input,p.Ti,RMSE_pd,RMSE_t
    ts_model = success_T(l,1);
    kappa_g_model = success_T(l,2);
    % cc_mask = ( seeds_stage2(:,1)==ts_model & seeds_stage2(:,2)==kappa_g_model );
    kappa_r_model=success_T(l,3);
    Rs_model=success_T(l,4);
    Rp_model=success_T(l,5);

    BSE_model = -4;
    % BSE_model = success_T(l,6); 
    Ubse_model = success_T(l,6); 
    Qc_model = success_T(l,7);
    dQc_model = success_T(l,8);
    Ti_model = success_T(l,9);
    eta_model = success_T(l,10);

    % Calculate the corresponding crustal growth pattern
%     [NGcc_model,Rcc_model,Gcc_model] = CC_growth_fun1(t,ts_model,tmax,...
%     Mcp,kappa_g_model,Rp_model,Rs_model,kappa_r_model);
    futures1(l)=parfeval(@CC_growth_fun1, 3, t,ts_model,tmax,...
     Mcp,kappa_g_model,Rp_model,Rs_model,kappa_r_model);

%     % save the results
%     NGcc_T(:,l) = NGcc_model;
%     Rcc_T(:,l) = Rcc_model;
%     Gcc_T(:,l) = Gcc_model;

    % Calculate the corresponding Temperature curve
%     [Tp_model,Tm_model,Qs_model,Qc_model,H_model,Ra_model,eta_out] = Thermal_model(t*1e9*yr_s,Ti_model,eta_model,Qc_model,dQc_model,...
%             BSE_model,Ubse_model,kappa_g_model,ts_model*1e9*yr_s,para_Scale,Hsf);
    futures2(l)=parfeval(@Thermal_model, 7, t*1e9*yr_s,Ti_model,eta_model,Qc_model,dQc_model,...
            BSE_model,Ubse_model,kappa_g_model,ts_model*1e9*yr_s,para_Scale,Hsf);

%     Tp_T(:,l) = Tp_model;
%     Tm_T(:,l) = Tm_model;
%     Qs_T(:,l) = Qs_model;
%     Qc_T(:,l) = Qc_model;
%     H_T(:,l) = H_model;
%     Ra_T(:,l) = Ra_model;

end % for l=1:itermax

elapsedTime = toc;  % s
disp(['used time = ', num2str(elapsedTime), ' s']);

for i = 1:size_T
    if mod(i,10)==0
        disp(['lT11=' num2str(i) ' of ' num2str(size_T)]);% keep track of the calculation
    end
    try        
        [NGcc_T(:,i),Rcc_T(:,i),Gcc_T(:,i)]=fetchOutputs(futures1(i)); 
        [Tp_T(:,i),Tm_T(:,i),Qs_T(:,i),Qc_T(:,i),H_T(:,i),Ra_T(:,i),~]=fetchOutputs(futures2(i));
    catch ME
        warning('Error in task %d: %s', i, ME.message);      
        NGcc_T(:,i)=nan(nt,1);
        Rcc_T(:,i)=nan(nt,1);
        Gcc_T(:,i)=nan(nt,1);
        Tp_T(:,i)=nan(nt,1);
        Tm_T(:,i)=nan(nt,1);
        Qs_T(:,i)=nan(nt,1);
        Qc_T(:,i)=nan(nt,1);
        H_T(:,i)=nan(nt,1);
        Ra_T(:,i)=nan(nt,1);

    end
end
elapsedTime = toc;  % s
disp(['used time = ', num2str(elapsedTime), ' s']);



% save
T_output.NG_cc = NGcc_T;
T_output.G_cc = Gcc_T;
T_output.R_cc = Rcc_T;
% T_output.F_T = F_T;
% T_output.S_T = S_T;
T_output.Tp_T=Tp_T;
T_output.Tm_T=Tm_T;
T_output.Qs_T=Qs_T;
T_output.Qc_T=Qc_T;
T_output.H_T=H_T;
T_output.Ra_T=Ra_T;

% save('thermal/thermal_output_for_figure.mat','T_output','-v7.3');






%% % AFTER Xe
% timestep too small, cannot run with parpool (storage limit)
% takes 230s
% 12.8G data, prefer generate than store
%% load observation 
% present day mantle processing rate
P_mor_pd = [1e14 10e14];%kg/yr
 
% Load in the observed water & Xe observation
Nom = 1;
Nom_range = 0.05;
water_obs1 = struct('NOM', Nom,'NOM_range',Nom_range);

% observation of Xe: present day
    Xe_obs = struct( 'Xe',       [4.3e5	9.2e5], ... % atoms/g
                     'Xe128v130',     [0.475	0.478], ...
                     'Xe128v132',     [0.069	0.071], ...
                     'Xe130v132',     [0.1445	0.1493], ...
                     'Xe131v132',     [0.7608	0.7786], ...
                     'Xe134v132',     [0.4082	0.4302], ...
                     'Xe136v132',     [0.3559	0.3835] ...
                    );

default_input_for_degassing;
t=Simulation_input.time_series;
t=t/yr_s/1e9;% to Ga
nt = length(t);
ntd=nt;
td=t;
%% full run

tic;
% size_Xe=12;

% futures1(size_Xe) = parallel.FevalFuture;

cc_Xe = nan(nt,size_Xe);
NGcc_Xe = nan(nt,size_Xe);
Rcc_Xe = nan(nt,size_Xe);
Gcc_Xe = nan(nt,size_Xe);

Tp_Xe = nan(ntd,size_Xe);
Tm_Xe = nan(ntd,size_Xe);
Qs_Xe = nan(ntd,size_Xe);
Qc_Xe = nan(ntd,size_Xe);
H_Xe = nan(ntd,size_Xe);
Ra_Xe = nan(ntd,size_Xe);
U_Xe = nan(ntd,size_Xe);

D0_cc=nan(ntd,size_Xe);
D0_mor=nan(ntd,size_Xe);
D0_p=nan(ntd,size_Xe);
P0_mor=nan(ntd,size_Xe);
P0_p=nan(ntd,size_Xe);

Xe = nan(ntd,size_Xe);
Xe128 = nan(ntd,size_Xe);
Xe131 = nan(ntd,size_Xe);
Xe132 = nan(ntd,size_Xe);
Xe134 = nan(ntd,size_Xe);
Xe136 = nan(ntd,size_Xe);

Pu = nan(ntd,size_Xe);
Ur = nan(ntd,size_Xe);



for l=1:size_Xe
    if mod(l,10)==0
        disp(['lXe0=' num2str(l) ' of ' num2str(size_Xe)]);% keep track of the calculation
    end
    
    % independent main variables
    param_success = struct ('time_cc_input', success_Xe(l,1)*1e9*yr_s,...
                            'kappa_gcc', success_Xe(l,2),...
                            'kappa_rcc', success_Xe(l,3),...
                            'Rs', success_Xe(l,4),...
                            'Rp', success_Xe(l,5),...
                            'BSE', -4,...% success_Xe(l,6),...
                            'Ubse', success_Xe(l,6),...% 1.273e-8,...
                            'Qc_pd_input', success_Xe(l,7),...
                            'dQc_input', success_Xe(l,8),...
                            'Ti', success_Xe(l,9),...
                            'eta_ref', success_Xe(l,10),...
                            'Xm_init', success_Xe(l,11),...
                            'Ms_init', success_Xe(l,12),...
                            'P_mor_pd', P_mor_pd,...
                            'Fr_w', success_Xe(l,13),...
                            'Xe_init', success_Xe(l,14),...
                            'Fr_Xe', success_Xe(l,15),...
                            'Fd_mor', success_Xe(l,16),...
                            'Fd_p', success_Xe(l,17),...
                            't',       t*1e9*yr_s, ... % change from Ga to s
                            'tmax',    tmax*1e9*yr_s, ... % change from Ga to s ====== above: time vector =======
                            'water_obs1', water_obs1,... % ====== above: watre model ====== 
                            'SwXe',  SwXe,...
                            'Xe_obs', Xe_obs ... % ====== above: Xe model ======
                             );

    [CC_growth_input,Thermal_input,H2O_input,Xe_input,Simulation_input,FIG,Fdmor_model,Fdp_model] = input_degassing(param_success);

    
        
        % Simulation_input.nsteps  = Simulation_input.end*[ 4601 ]; % number of timesteps
        output_model = Degassing_model(Simulation_input,...
                     Thermal_input,CC_growth_input,H2O_input,Xe_input,Fdmor_model,Fdp_model);
%         futures1(l)=parfeval(@Degassing_model, 1, Simulation_input,...
%                       Thermal_input,CC_growth_input,H2O_input,Xe_input,Fdmor_model,Fdp_model);

%         Rcc_model=output_model.R_cc*1e9*yr_s;
%         NGcc1_model=output_model.NG_cc*1e9*yr_s;
%         Gcc_model=Rcc_model+NGcc1_model;
%         NGcc_model=output_model.cc*Mcp;
%         ts_model=success_Xe(l,1);
%         kappa_g_model=success_Xe(l,2);
%         Rp_model=success_Xe(l,3);
%         Rs_model=success_Xe(l,4);
%         kappa_r_model=success_Xe(l,5);
%         frw_model=Krw_Xe(l);
% 
%         [NGcc_model,Rcc_model,Gcc_model,Krw_model] = CC_growth_fun2(t,ts_model,tmax,...
%         Mcp,kappa_g_model,Rp_model,Rs_model,kappa_r_model,frw_model);
%         [F_model, S_model, m_tp,m,Krw_real_model,s] = Formation_surface_age_fun(...
%             t, Gcc_model, Rcc_model, NGcc_model, Krw_model);


        
        % save the results
        cc_Xe(:,l) = output_model.cc;
        NGcc_Xe(:,l) = output_model.NG_cc;
        Rcc_Xe(:,l) = output_model.R_cc;
%         F_Xe(:,l) = F_model;
%         S_Xe(:,l) = S_model;

        Tm_Xe(:,l) = output_model.T;
        Qs_Xe(:,l) = output_model.Qs;
        Qc_Xe(:,l) = output_model.Qc;
        H_Xe(:,l) = output_model.H;
        Ra_Xe(:,l) = output_model.Ra;
        U_Xe(:,l) = output_model.U; % velocity

        D0_mor(:,l)=output_model.D0_mor;
        D0_cc(:,l)=output_model.D0_cc;
        D0_p(:,l)=output_model.D0_p;
        P0_mor(:,l)=output_model.D0_mor./param_success.Fd_mor;
        P0_p(:,l)=output_model.D0_p./param_success.Fd_p;

        Xe(:,l) =  output_model.Xe;
        Xe128(:,l) = output_model.Xe_Atm;
        Xe131(:,l) = output_model.Xe131;
        Xe132(:,l) = output_model.Xe132;
        Xe134(:,l) = output_model.Xe134;
        Xe136(:,l) = output_model.Xe136;
        Pu(:,l) = output_model.Put;
        Ur(:,l) = output_model.Urt;

        % problem: Formation_surface_age_fun only for uniform timestep
%       [NGcc_model,Rcc_model,Gcc_model,Krw_model] = CC_growth_fun2(t,ts_model,tmax,...
%          Mcp,kappa_g_model,Rp_model,Rs_model,kappa_r_model,frw_model);
%       [F_model, S_model, m_tp,m,Krw_real_model,s] = Formation_surface_age_fun(...
%           param_success.t, Mud_model, Mdd_model, Mc_model, Krw_model);
% 



end % for l=1:itermax
% elapsedTime = toc;  % s
% disp(['used time = ', num2str(elapsedTime), ' s']);

elapsedTime = toc;  % s
disp(['used time = ', num2str(elapsedTime), ' s']);

% %%
% dt = 0.001;% length of each timestep
% t = 0:dt:tmax;
% nt = length(t);% number of timesteps
% t = t';
% 
% F_Xe = nan(nt,size_Xe);
% S_Xe = nan(nt,size_Xe);
% 
% tic;
% size_Xe=1;
% for l=1:size_Xe
%     if mod(l,10)==0
%         disp(['lXe0=' num2str(l) ' of ' num2str(size_Xe)]);% keep track of the calculation
%     end
%        
%         ts_model=success_Xe(l,1);
%         kappa_g_model=success_Xe(l,2);
%         Rp_model=success_Xe(l,3);
%         Rs_model=success_Xe(l,4);
%         kappa_r_model=success_Xe(l,5);
%         frw_model=Krw_Xe(l);
% 
%         [NGcc_model,Rcc_model,Gcc_model,Krw_model] = CC_growth_fun2(t,ts_model,tmax,...
%         Mcp,kappa_g_model,Rp_model,Rs_model,kappa_r_model,frw_model);
%         [F_model, S_model, m_tp,m,Krw_real_model,s] = Formation_surface_age_fun(...
%             t, Gcc_model, Rcc_model, NGcc_model, Krw_model);
%        
%         % save the results     
%         F_Xe(:,l) = F_model;
%         S_Xe(:,l) = S_model;
% 
% end % for l=1:itermax
% % elapsedTime = toc;  % s
% % disp(['used time = ', num2str(elapsedTime), ' s']);

% elapsedTime = toc;  % s
% disp(['used time = ', num2str(elapsedTime), ' s']);

% for i = 1:size_Xe
%     if mod(i,10)==0
%         disp(['lXe1=' num2str(i) ' of ' num2str(size_Xe)]);% keep track of the calculation
%     end
%     try
%         % [Tp_T(:,i),Tm_T(:,i),Qs_T(:,i),Qc_T(:,i),H_T(:,i),Ra_T(:,i),~]=fetchOutputs(futures2(i));
%         output_model=fetchOutputs(futures1(i));
%         cc_Xe(:,i) = output_model.cc;
%         NGcc_Xe(:,i) = output_model.NG_cc;
%         Rcc_Xe(:,i) = output_model.R_cc;
%         Tm_Xe(:,i) = output_model.T;
%         Qs_Xe(:,i) = output_model.Qs;
%         Qc_Xe(:,i) = output_model.Qc;
%         H_Xe(:,i) = output_model.H;
%         Ra_Xe(:,i) = output_model.Ra;
%         U_Xe(:,i) = output_model.U; % velocity
% 
%         Xe(:,i) =  output_model.Xe;
%         Xe128(:,i) = output_model.Xe_Atm;
%         Xe131(:,i) = output_model.Xe131;
%         Xe132(:,i) = output_model.Xe132;
%         Xe134(:,i) = output_model.Xe134;
%         Xe136(:,i) = output_model.Xe136;
%         Pu(:,i) = output_model.Put;
%         Ur(:,i) = output_model.Urt;
% 
%     catch ME
%         warning('Error in task %d: %s', i, ME.message);      
%         cc_Xe(:,i) = nan(nt,1);
%         NGcc_Xe(:,i) = nan(nt,1);
%         Rcc_Xe(:,i) = nan(nt,1);
%         Tm_Xe(:,i) = nan(nt,1);
%         Qs_Xe(:,i) = nan(nt,1);
%         Qc_Xe(:,i) = nan(nt,1);
%         H_Xe(:,i) = nan(nt,1);
%         Ra_Xe(:,i) = nan(nt,1);
%         U_Xe(:,i) =nan(nt,1); % velocity
% 
%         Xe(:,i) =  nan(nt,1);
%         Xe128(:,i) = nan(nt,1);
%         Xe131(:,i) = nan(nt,1);
%         Xe132(:,i) = nan(nt,1);
%         Xe134(:,i) = nan(nt,1);
%         Xe136(:,i) = nan(nt,1);
%         Pu(:,i) = nan(nt,1);
%         Ur(:,i) = nan(nt,1);
%     end
% end
% elapsedTime = toc;  % s
% disp(['used time = ', num2str(elapsedTime), ' s']);


% fr kg/s to kg/Ga
NGcc_Xe = NGcc_Xe*1e9*yr_s;
Rcc_Xe = Rcc_Xe*1e9*yr_s;
Gcc_Xe = NGcc_Xe+Rcc_Xe;
% Tm to Tp
Tp_Xe = Tm_Xe*exp(-alpha_ra*g*1445e3/Cp)-273; %C
%



% save
Xe_output.cc_Xe=cc_Xe  ;
Xe_output.NG_cc=NGcc_Xe  ;
Xe_output.G_cc=Gcc_Xe  ;
Xe_output.R_cc=Rcc_Xe  ;
% Xe_output.F_Xe=F_Xe  ;
% Xe_output.S_Xe=S_Xe  ;
Xe_output.Tm_Xe=Tm_Xe;
Xe_output.Qs_Xe=Qs_Xe;
Xe_output.Qc_Xe=Qc_Xe;
Xe_output.H_Xe=H_Xe;
Xe_output.Ra_Xe=Ra_Xe;
Xe_output.U_Xe=U_Xe ; % velocity

Xe_output.D0_p=D0_p;
Xe_output.D0_mor=D0_mor;
Xe_output.D0_cc=D0_cc;

Xe_output.Xe=Xe;
Xe_output.Xe128=Xe128;
Xe_output.Xe131=Xe131;
Xe_output.Xe132=Xe132;
Xe_output.Xe134=Xe134;
Xe_output.Xe136=Xe136;
Xe_output.Pu=Pu;
Xe_output.Ur=Ur;

% save('Xe/Xe_output_for_figure.mat','Xe_output','-v7.3');
% save('Xe/Tp_Xe','Tp_Xe','-v7.3');


%% ex run

tic;
l=randi(size_Xe);
size_Xe=5;

kappa_gcc_ex=[success_Xe(l,2),-1,success_Xe(l,2),success_Xe(l,2),success_Xe(l,2)];
kappa_rcc_ex=[success_Xe(l,3),success_Xe(l,3),-3,success_Xe(l,3),success_Xe(l,3)];
Rs_ex=[success_Xe(l,4),success_Xe(l,4),success_Xe(l,4),0,success_Xe(l,3)];
Ubse_ex=[success_Xe(l,6),success_Xe(l,6),success_Xe(l,6),success_Xe(l,6),20e-9];

% futures1(size_Xe) = parallel.FevalFuture;


cc_Xe = nan(nt,size_Xe);
NGcc_Xe = nan(nt,size_Xe);
Rcc_Xe = nan(nt,size_Xe);
Gcc_Xe = nan(nt,size_Xe);

Tp_Xe = nan(ntd,size_Xe);
Tm_Xe = nan(ntd,size_Xe);
Qs_Xe = nan(ntd,size_Xe);
Qc_Xe = nan(ntd,size_Xe);
H_Xe = nan(ntd,size_Xe);
Ra_Xe = nan(ntd,size_Xe);
U_Xe = nan(ntd,size_Xe);

D0_mor=nan(ntd,size_Xe);
D0_cc=nan(ntd,size_Xe);
D0_p=nan(ntd,size_Xe);

Xe = nan(ntd,size_Xe);
Xe128 = nan(ntd,size_Xe);
Xe131 = nan(ntd,size_Xe);
Xe132 = nan(ntd,size_Xe);
Xe134 = nan(ntd,size_Xe);
Xe136 = nan(ntd,size_Xe);

Pu = nan(ntd,size_Xe);
Ur = nan(ntd,size_Xe);



for k=1:size_Xe
    disp(['lXe0=' num2str(l) ' of ' num2str(size_Xe)]);% keep track of the calculation
    
    
    % independent main variables
    param_success = struct ('time_cc_input', success_Xe(l,1)*1e9*yr_s,...
                            'kappa_gcc', kappa_gcc_ex(k),...% success_Xe(l,2),...
                            'kappa_rcc', kappa_rcc_ex(k),...
                            'Rs', Rs_ex(k),...
                            'Rp', success_Xe(l,5),...
                            'BSE', -4,...% success_Xe(l,6),...
                            'Ubse', Ubse_ex(k),...% success_Xe(l,6),...% 1.273e-8,...
                            'Qc_pd_input', success_Xe(l,7),...
                            'dQc_input', success_Xe(l,8),...
                            'Ti', success_Xe(l,9),...
                            'eta_ref', success_Xe(l,10),...
                            'Xm_init', success_Xe(l,11),...
                            'Ms_init', success_Xe(l,12),...
                            'P_mor_pd', P_mor_pd,...
                            'Fr_w', success_Xe(l,13),...
                            'Xe_init', success_Xe(l,14),...
                            'Fr_Xe', success_Xe(l,15),...
                            'Fd_mor', success_Xe(l,16),...
                            'Fd_p', success_Xe(l,17),...
                            't',       t*1e9*yr_s, ... % change from Ga to s
                            'tmax',    tmax*1e9*yr_s, ... % change from Ga to s ====== above: time vector =======
                            'water_obs1', water_obs1,... % ====== above: watre model ====== 
                            'SwXe',  SwXe,...
                            'Xe_obs', Xe_obs ... % ====== above: Xe model ======
                             );

    [CC_growth_input,Thermal_input,H2O_input,Xe_input,Simulation_input,FIG,Fdmor_model,Fdp_model] = input_degassing(param_success);

    
        
        % Simulation_input.nsteps  = Simulation_input.end*[ 4601 ]; % number of timesteps
        output_model = Degassing_model(Simulation_input,...
                     Thermal_input,CC_growth_input,H2O_input,Xe_input,Fdmor_model,Fdp_model);
        
        % save the results
        cc_Xe(:,k) = output_model.cc;
        NGcc_Xe(:,k) = output_model.NG_cc;
        Rcc_Xe(:,k) = output_model.R_cc;
%         F_Xe(:,k) = F_model;
%         S_Xe(:,k) = S_model;

        Tm_Xe(:,k) = output_model.T;
        Qs_Xe(:,k) = output_model.Qs;
        Qc_Xe(:,k) = output_model.Qc;
        H_Xe(:,k) = output_model.H;
        Ra_Xe(:,k) = output_model.Ra;
        U_Xe(:,k) = output_model.U; % velocity

        D0_mor(:,k)=output_model.D0_mor;
        D0_cc(:,k)=output_model.D0_cc;
        D0_p(:,k)=output_model.D0_p;

        Xe(:,k) =  output_model.Xe;
        Xe128(:,k) = output_model.Xe_Atm;
        Xe131(:,k) = output_model.Xe131;
        Xe132(:,k) = output_model.Xe132;
        Xe134(:,k) = output_model.Xe134;
        Xe136(:,k) = output_model.Xe136;
        Pu(:,k) = output_model.Put;
        Ur(:,k) = output_model.Urt;



end 

elapsedTime = toc;  % s
disp(['used time = ', num2str(elapsedTime), ' s']);


% fr kg/s to kg/Ga
NGcc_Xe = NGcc_Xe*1e9*yr_s;
Rcc_Xe = Rcc_Xe*1e9*yr_s;
Gcc_Xe = NGcc_Xe+Rcc_Xe;
% Tm to Tp
Tp_Xe = Tm_Xe*exp(-alpha_ra*g*1445e3/Cp)-273; %C

% save
Xe_output.cc_Xe=cc_Xe  ;
Xe_output.NG_cc=NGcc_Xe  ;
Xe_output.G_cc=Gcc_Xe  ;
Xe_output.R_cc=Rcc_Xe  ;
% Xe_output.F_Xe=F_Xe  ;
% Xe_output.S_Xe=S_Xe  ;
Xe_output.Tm_Xe=Tm_Xe;
Xe_output.Qs_Xe=Qs_Xe;
Xe_output.Qc_Xe=Qc_Xe;
Xe_output.H_Xe=H_Xe;
Xe_output.Ra_Xe=Ra_Xe;
Xe_output.U_Xe=U_Xe ; % velocity

Xe_output.D0_p=D0_p;
Xe_output.D0_mor=D0_mor;
Xe_output.D0_cc=D0_cc;

Xe_output.Xe=Xe;
Xe_output.Xe128=Xe128;
Xe_output.Xe131=Xe131;
Xe_output.Xe132=Xe132;
Xe_output.Xe134=Xe134;
Xe_output.Xe136=Xe136;
Xe_output.Pu=Pu;
Xe_output.Ur=Ur;

% save('Xe/Xe_output_for_figure_ex.mat','Xe_output','-v7.3');
% save('Xe/Tp_Xe','Tp_Xe','-v7.3');

