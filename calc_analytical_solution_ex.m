%% % sup figures


%% % analytical solution vs numerical solution

%% load specified numerical results
Xe_output_struct=load('Xe/out_Xe_lhs_new_3e7_test_fixK_U_MC_lowerT_alot_Rac1100_new_lowf.mat');
results=Xe_output_struct.results;
col_names=Xe_output_struct.col_names;

crit_misfit_Pmor = 1;
crit_misfit_w = 1;
crit_misfit_xec = 1;
crit_misfit_xe128v130 = 1; 
crit_misfit_xe128v132 = 1; 
crit_misfit_xe130v132 = 1; 
crit_misfit_xe131v132 = 1; 
crit_misfit_xe134v132 = 1; 
crit_misfit_xe136v132 = 1; 

U_idx=find(strcmp(col_names, 'Ubse'));
ErPmor_idx = find(strcmp(col_names, 'Er_Pmor_pd'));
rmsew_idx = find(strcmp(col_names, 'RMSE_water'));
Erxec_idx = find(strcmp(col_names, 'Er_xec'));
Erxe128v130_idx = find(strcmp(col_names, 'Er_xe128v130'));
Erxe128v132_idx = find(strcmp(col_names, 'Er_xe128v132'));
Erxe130v132_idx = find(strcmp(col_names, 'Er_xe130v132'));
Erxe131v132_idx = find(strcmp(col_names, 'Er_xe131v132'));
Erxe134v132_idx = find(strcmp(col_names, 'Er_xe134v132'));
Erxe136v132_idx = find(strcmp(col_names, 'Er_xe136v132'));

success_mask = (    abs( results(:, Erxec_idx) ) <= crit_misfit_xec & ...
                    abs( results(:, Erxe128v130_idx) ) <= crit_misfit_xe128v130 & ... % results(:, Erxe128v132_idx) <= crit_misfit_xe128v132 & ...                    
                    results(:, U_idx)==20e-9);
%                     abs( results(:, Erxe130v132_idx) ) <= crit_misfit_xe130v132 & ...
%                     abs( results(:, Erxe131v132_idx) ) <= crit_misfit_xe131v132 & ...
%                     abs( results(:, Erxe134v132_idx) ) <= crit_misfit_xe134v132 & ...
%                     abs( results(:, Erxe136v132_idx) ) <= crit_misfit_xe136v132 & ...
%                     abs( results(:, ErPmor_idx) ) <= crit_misfit_Pmor ); %& ...
                    %results(:, rmsew_idx) <= crit_misfit_w & ...
                 
Xe_ex_U=results(success_mask,:);

%% load mean numerical solution

Xe_output_margin_struct=load('Xe/Xe_output_for_figure_margin_new_lowf.mat');
Xe_output_margin=Xe_output_margin_struct.Xe_output_margin;
Xe_output_ex_struct=load('Xe/Xe_output_for_figure_ex_new_lowf.mat');
Xe_output_ex=Xe_output_ex_struct.Xe_output;

Xe_50=Xe_output_margin.Xe_50;
Xe128_50=Xe_output_margin.Xe128_50;
Xe131_50=Xe_output_margin.Xe131_50;
Xe132_50=Xe_output_margin.Xe132_50;
Xe134_50=Xe_output_margin.Xe134_50;
Xe136_50=Xe_output_margin.Xe136_50;

Xe128r130_50=Xe_output_margin.Xe128r130_50;

Dcc_50=Xe_output_margin.D0_cc_50;
Dm_50=Xe_output_margin.D0_mor_50;
Dp_50=Xe_output_margin.D0_p_50;

Xe_ex=Xe_output_ex.Xe;
Dcc_ex=Xe_output_ex.D0_cc;
Dm_ex=Xe_output_ex.D0_mor;
Dp_ex=Xe_output_ex.D0_p;


%% load initial controlling params in the analytical solution
par_solver_new;
default_input_for_degassing;
t_n=Simulation_input.time_series;
t_n=t_n/yr_s/1e9;
Xe_obs = struct( 'Xe',       [4.3e5	9.2e5], ... % atoms/g
                 'Xe128v130',     [0.475	0.478], ...
                 'Xe128v132',     [0.069	0.071], ...
                 'Xe130v132',     [0.1445	0.1493], ...
                 'Xe131v132',     [0.7608	0.7786], ...
                 'Xe134v132',     [0.4082	0.4302], ...
                 'Xe136v132',     [0.3559	0.3835] ...
                );

Xe_obs1 = struct( 'Xe',       [4.3e5	9.2e5], ... % atoms/g
                 'Xe128v130',     [0.475	0.478], ...
                 'Xe128v132',     [0.069	0.071], ...
                 'Xe132v130',     [1/0.1493 1/0.1445], ...
                 'Xe131v130',     [0.7608/0.1493	0.7786/0.1445], ...
                 'Xe134v130',     [0.4082/0.1493	0.4302/0.1445], ...
                 'Xe136v130',     [0.3559/0.1493	0.3835/0.1445] ...
                );


% Xe in regassing (in hydrous mineral) atmos/g 8e6 used in this 1.262e-9 (=5.85e6atoms/g)
Xe_slab=Xe_input.XesAtm_init; 

lamU=Xe_input.lamUr; % [yr-1]
lamPu=Xe_input.lamPu; % [yr-1]
% YPu=7e-5;%136Xe
% YUr=3.43e-8;%136Xe
% YPu132=YPu/1.738;%132Xe
% YUr132=YUr/1.120;%132Xe
% YPu131=YPu132*0.1449;%131Xe
% YUr131=YUr132*0.2777;%131Xe
% YPu134=YPu132*1.437;%134Xe
% YUr134=YUr132*1.041;%134Xe
% 

M=M;
tau_U=1/lamU/1e9; % [yr]->[Ga]
tau_Pu=1/lamPu/1e9; % [yr]->[Ga]

% give solution
% initial: en; atm: pd
% m128=(r128i-r128pd)/(r128i-r128a);
m128_range_en=[(Xeren-Xe_obs.Xe128v130(2))/(Xeren-Xerpa),(Xeren-Xe_obs.Xe128v130(1))/(Xeren-Xerpa)];
% Xe_R=Xe_pd*m128
Xe_R_range_en=[Xe_obs.Xe(1)*m128_range_en(1),Xe_obs.Xe(2)*m128_range_en(2)]; % [atmos/g]
% tau_D=t_pd/(ln(Xe_i/(Xe_pd-Xe_R)))=t_pd/(ln(Xe_i/(Xe_pd-Xe_pd*m128))=t_pd/(ln(Xe_i/(Xe_pd*(1-m128)))
% tau_D=t_pd/(ln(Xe_i/Xe_pd/(1-m128)))
%tau_D_range_en=4.6./[log(3.2e8/Xe_obs.Xe(1)/(1-m128_range_en(2))),log(3.2e7/Xe_obs.Xe(2)/(1-m128_range_en(1)))]; % [Gyr]
tau_D_range_en=4.6./[log(3.2e7/Xe_obs.Xe(1)/(1-m128_range_en(2))),log(3.2e7/Xe_obs.Xe(2)/(1-m128_range_en(1)))]; % [Gyr]

m128_range=m128_range_en;
Xe_R_range=Xe_R_range_en;
tau_D_range=tau_D_range_en;

% give solution
% initial: en; atm: pd
% m128=(r128i-r128pd)/(r128i-r128a);
m128_range_cc=[(Xercc-Xe_obs.Xe128v130(2))/(Xercc-Xerpa),(Xercc-Xe_obs.Xe128v130(1))/(Xercc-Xerpa)];
% Xe_R=Xe_pd*m128
Xe_R_range_cc=[Xe_obs.Xe(1)*m128_range_cc(1),Xe_obs.Xe(2)*m128_range_cc(2)]; % [atmos/g]
% tau_D=t_pd/(ln(Xe_i/(Xe_pd-Xe_R)))=t_pd/(ln(Xe_i/(Xe_pd-Xe_pd*m128))=t_pd/(ln(Xe_i/(Xe_pd*(1-m128)))
% tau_D=t_pd/(ln(Xe_i/Xe_pd/(1-m128)))
%tau_D_range_cc=4.6./[log(3.2e8/Xe_obs.Xe(1)/(1-m128_range_cc(2))),log(3.2e7/Xe_obs.Xe(2)/(1-m128_range_cc(1)))]; % [Gyr]
tau_D_range_cc=4.6./[log(3.2e7/Xe_obs.Xe(1)/(1-m128_range_cc(2))),log(3.2e7/Xe_obs.Xe(2)/(1-m128_range_cc(1)))]; % [Gyr]


% U_i=Xe_pd*(rf_pd-rf_a*m128-rf_i*(1-m128))/(-KPu*(Pu_i/U_i)*exp(-t_pd/tau_D)+KU*exp(-t_pd/tau_U));
% U_i=Xe_pd*(rf_pd-rf_a*m128-rf_i*(1-m128))/(-KPu*(Pu_i/U_i)*(Xe_pd/Xe_i*(1-m128))+KU*exp(-t_pd/tau_U))
R0Pu_Ur=0.0068;

%%
m128_step=(m128_range(2)-m128_range(1))/(1e2-1);
m128_test=m128_range(1):m128_step:m128_range(2);
tau_D_min_test=4.6./log(3.2e7/Xe_obs.Xe(1)./(1-m128_test)); % [Gyr] min Xei & min Xepd
tau_D_max_test=4.6./log(3.2e8/Xe_obs.Xe(2)./(1-m128_test)); % [Gyr] % max Xei & max Xepd
% tau_D_max_test=4.6./log(3.2e8/Xe_obs.Xe(2)./(1-m128_test)); % [Gyr] % max Xei & max Xepd
% tau_D_min_test=min(4.6./log(1./(1-m128_test))) % [Gyr] 
% tau_D_max_test=max(4.6./log(1./(1-m128_test))) % [Gyr] % 


KPu131_min_test=YPu131./(tau_Pu./tau_D_min_test-1);
KU131_min_test=YUr131./(tau_U./tau_D_min_test-1);
KPu131_max_test=YPu131./(tau_Pu./tau_D_max_test-1);
KU131_max_test=YUr131./(tau_U./tau_D_max_test-1);


KPu=KPu131_min_test;KU=KU131_min_test;
% [Ubse_min_131,loc1]=min( (Xe_obs.Xe(1).*(Xe_obs1.Xe131v130(1)-Xe131rpa.*m128_test-Xe131ren.*(1-m128_test))...
%     ./(-KPu.*R0Pu_Ur.*exp(-4.6./tau_D_min_test)+KU.*exp(-4.6/tau_U))) ...
%     *238/NA*1e9*exp(-Xe_input.lamUr*t_pd/yr_s) );
[Ubse_min_131,loc1]=min( (Xe_obs.Xe(1).*(Xe_obs1.Xe131v130(1)-Xe131rpa.*m128_test-Xe131ren.*(1-m128_test))...
    ./(-KPu.*R0Pu_Ur.*exp(-4.6./tau_D_min_test)+KU.*(exp(-4.6/tau_U)-exp(-4.6./tau_D_min_test) ) ) ...
    *238/NA*1e9*exp(-Xe_input.lamUr*t_pd/yr_s) ) )
KPu=KPu131_max_test;KU=KU131_max_test;
% [Ubse_max_131,loc2]=max( (Xe_obs.Xe(2).*(Xe_obs1.Xe131v130(2)-Xe131rpa.*m128_test-Xe131ren.*(1-m128_test))...
%     ./(-KPu.*R0Pu_Ur.*exp(-4.6./tau_D_max_test)+KU.*exp(-4.6/tau_U)))...
%     *238/NA*1e9*exp(-Xe_input.lamUr*t_pd/yr_s) )
[Ubse_max_131,loc2]=max( (Xe_obs.Xe(2).*(Xe_obs1.Xe131v130(2)-Xe131rpa.*m128_test-Xe131ren.*(1-m128_test))...
    ./(-KPu.*R0Pu_Ur.*exp(-4.6./tau_D_max_test)+KU.*(exp(-4.6/tau_U)- exp(-4.6./tau_D_max_test)) ))...
    *238/NA*1e9*exp(-Xe_input.lamUr*t_pd/yr_s) )

KPu132_min_test=YPu132./(tau_Pu./tau_D_min_test-1);
KU132_min_test=YUr132./(tau_U./tau_D_min_test-1);
KPu132_max_test=YPu132./(tau_Pu./tau_D_max_test-1);
KU132_max_test=YUr132./(tau_U./tau_D_max_test-1);

KPu=KPu132_min_test;KU=KU132_min_test;
[Ubse_min_132,loc1]=min( (Xe_obs.Xe(1).*(Xe_obs1.Xe132v130(1)-Xe132rpa.*m128_test-Xe132ren.*(1-m128_test))...
    ./(-KPu.*R0Pu_Ur.*exp(-4.6./tau_D_min_test)+KU.*exp(-4.6/tau_U))) ...
    *238/NA*1e9*exp(-Xe_input.lamUr*t_pd/yr_s) )  
KPu=KPu132_max_test;KU=KU132_max_test;
[Ubse_max_132,loc2]=max( (Xe_obs.Xe(2).*(Xe_obs1.Xe132v130(2)-Xe132rpa.*m128_test-Xe132ren.*(1-m128_test))...
    ./(-KPu.*R0Pu_Ur.*exp(-4.6./tau_D_max_test)+KU.*exp(-4.6/tau_U)))...
    *238/NA*1e9*exp(-Xe_input.lamUr*t_pd/yr_s) )

KPu134_min_test=YPu134./(tau_Pu./tau_D_min_test-1);
KU134_min_test=YUr134./(tau_U./tau_D_min_test-1);
KPu134_max_test=YPu134./(tau_Pu./tau_D_max_test-1);
KU134_max_test=YUr134./(tau_U./tau_D_max_test-1);

KPu=KPu134_min_test;KU=KU134_min_test;
[Ubse_min_134,loc1]=min( (Xe_obs.Xe(1).*(Xe_obs1.Xe134v130(1)-Xe134rpa.*m128_test-Xe134ren.*(1-m128_test))...
    ./(-KPu.*R0Pu_Ur.*exp(-4.6./tau_D_min_test)+KU.*exp(-4.6/tau_U))) ...
    *238/NA*1e9*exp(-Xe_input.lamUr*t_pd/yr_s) )  
KPu=KPu134_max_test;KU=KU134_max_test;
[Ubse_max_134,loc2]=max( (Xe_obs.Xe(2).*(Xe_obs1.Xe134v130(2)-Xe134rpa.*m128_test-Xe134ren.*(1-m128_test))...
    ./(-KPu.*R0Pu_Ur.*exp(-4.6./tau_D_max_test)+KU.*exp(-4.6/tau_U)))...
    *238/NA*1e9*exp(-Xe_input.lamUr*t_pd/yr_s) )

KPu136_min_test=YPu./(tau_Pu./tau_D_min_test-1);
KU136_min_test=YUr./(tau_U./tau_D_min_test-1);
KPu136_max_test=YPu./(tau_Pu./tau_D_max_test-1);
KU136_max_test=YUr./(tau_U./tau_D_max_test-1);

KPu=KPu136_min_test;KU=KU136_min_test;
[Ubse_min_136,loc1]=min( (Xe_obs.Xe(1).*(Xe_obs1.Xe136v130(1)-Xe136rpa.*m128_test-Xe136ren.*(1-m128_test))...
    ./(-KPu.*R0Pu_Ur.*exp(-4.6./tau_D_min_test)+KU.*exp(-4.6/tau_U))) ...
    *238/NA*1e9*exp(-Xe_input.lamUr*t_pd/yr_s) )  
KPu=KPu136_max_test;KU=KU136_max_test;
[Ubse_max_136,loc2]=max( (Xe_obs.Xe(2).*(Xe_obs1.Xe136v130(2)-Xe136rpa.*m128_test-Xe136ren.*(1-m128_test))...
    ./(-KPu.*R0Pu_Ur.*exp(-4.6./tau_D_max_test)+KU.*exp(-4.6/tau_U)))...
    *238/NA*1e9*exp(-Xe_input.lamUr*t_pd/yr_s) )


% KPu=KPu131;KU=KU131;
% Ubse_range131=[0,...
%     Xe_obs.Xe(2)*(Xe_obs1.Xe131v130(2)-Xe131rpa*m128_range(1)-Xe131ren*(1-m128_range(1)))...
%     /(-KPu*R0Pu_Ur*Xe_obs.Xe(2)/3.2e8*(1-m128_range(1))+KU*exp(-4.6/tau_U))]...
%     *238/NA*1e9*exp(-Xe_input.lamUr*t_pd/yr_s); % [atmos/g] -> [ppb]
% KPu=KPu132;KU=KU132;
% Ubse_range132=[Xe_obs.Xe(1)*(Xe_obs1.Xe132v130(1)-Xe132rpa*m128_range(2)-Xe132ren*(1-m128_range(2)))...
%     /(-KPu*R0Pu_Ur*Xe_obs.Xe(1)/3.2e7*(1-m128_range(2))+KU*exp(-4.6/tau_U)),...
%     Xe_obs.Xe(2)*(Xe_obs1.Xe132v130(2)-Xe132rpa*m128_range(1)-Xe132ren*(1-m128_range(1)))...
%     /(-KPu*R0Pu_Ur*Xe_obs.Xe(2)/3.2e7*(1-m128_range(1))+KU*exp(-4.6/tau_U))]...
%     *238/NA*1e9*exp(-Xe_input.lamUr*t_pd/yr_s); % [atmos/g] -> [ppb]
% KPu=KPu134;KU=KU134;
% Ubse_range134=[Xe_obs.Xe(1)*(Xe_obs1.Xe134v130(1)-Xe134rpa*m128_range(1)-Xe134ren*(1-m128_range(1)))...
%     /(-KPu*R0Pu_Ur*Xe_obs.Xe(1)/3.2e7*(1-m128_range(1))+KU*exp(-4.6/tau_U)),...
%     Xe_obs.Xe(2)*(Xe_obs1.Xe134v130(2)-Xe134rpa*m128_range(2)-Xe134ren*(1-m128_range(2)))...
%     /(-KPu*R0Pu_Ur*Xe_obs.Xe(2)/3.2e8*(1-m128_range(2))+KU*exp(-4.6/tau_U))]...
%     *238/NA*1e9*exp(-Xe_input.lamUr*t_pd/yr_s); % [atmos/g] -> [ppb]
% KPu=KPu136;KU=KU136;
% Ubse_range136=[Xe_obs.Xe(1)*(Xe_obs1.Xe136v130(1)-Xe136rpa*m128_range(1)-Xe136ren*(1-m128_range(1)))...
%     /(-KPu*R0Pu_Ur*Xe_obs.Xe(1)/3.2e7*(1-m128_range(1))+KU*exp(-4.6/tau_U)),...
%     Xe_obs.Xe(2)*(Xe_obs1.Xe136v130(2)-Xe136rpa*m128_range(2)-Xe136ren*(1-m128_range(2)))...
%     /(-KPu*R0Pu_Ur*Xe_obs.Xe(2)/3.2e8*(1-m128_range(2))+KU*exp(-4.6/tau_U))]...
%     *238/NA*1e9*exp(-Xe_input.lamUr*t_pd/yr_s); % [atmos/g] -> [ppb]

%%
% give solution result
% D= % 4~8e24 kg/Ga
% R= % (0.8~0.9)*Xe_pd/Xe_slab*D
% for 128Xe/130Xe, larger D should refer to smaller R
D=5e24; % 5.128
R=0.85*Xe_obs.Xe(2)/Xe_slab*D;
%
% U=
U238bse=13e-9; % [ppb]->[g/g]
U_i=U238bse*exp(Xe_input.lamUr*t_pd/yr_s)/238*NA; % [g/g]->[atmos/g]
Pu_i=U_i*Xe_input.R0Pu_Ur; % [atmos/g]


% calculate dependent params
Xe130_R=Xe_slab*R/D;
% pd atm
Xe128_R1=Xe130_R*Xerpa;
Xe131_R1=Xe130_R*Xe131rpa;
Xe132_R1=Xe130_R*Xe132rpa;
Xe134_R1=Xe130_R*Xe134rpa;
Xe136_R1=Xe130_R*Xe136rpa;
% initial atm
Xe128_R2=Xe130_R*Xe_input.Xes_init;
Xe131_R2=Xe130_R*Xe_input.Xe131r130atm_init;
Xe132_R2=Xe130_R*Xe_input.Xe132r130atm_init;
Xe134_R2=Xe130_R*Xe_input.Xe134r130atm_init;
Xe136_R2=Xe130_R*Xe_input.Xe136r130atm_init;

Xe130_i=Xe_input.Xe_init;
Xe128_i=Xe_input.XeAtm_init;
Xe131_i=Xe_input.Xe131_init;
Xe132_i=Xe_input.Xe132_init;
Xe134_i=Xe_input.Xe134_init;
Xe136_i=Xe_input.Xe136_init;


tau_D=M/D; % [Ga]

KPu130=0;
KU130=0;
KPu128=0;
KU128=0;
KPu131=YPu131/(tau_Pu/tau_D-1);
KU131=YUr131/(tau_U/tau_D-1);
KPu132=YPu132/(tau_Pu/tau_D-1);
KU132=YUr132/(tau_U/tau_D-1);
KPu134=YPu134/(tau_Pu/tau_D-1);
KU134=YUr134/(tau_U/tau_D-1);
KPu136=YPu/(tau_Pu/tau_D-1);
KU136=YUr/(tau_U/tau_D-1);



%% % analytical solution
t=0:0.01:4.6;

Xe_R=Xe130_R;Xe_i=Xe130_i;
KPu=KPu130;KU=KU130;
Xe130=(Xe_i-Xe_R-KPu*Pu_i-KU*U_i)*exp(-t./tau_D)+Xe_R...
        +KPu*Pu_i*exp(-t./tau_Pu)+KU*U_i*exp(-t./tau_U);

Xe_R=Xe128_R1;Xe_i=Xe128_i;
KPu=KPu128;KU=KU128;
Xe128=(Xe_i-Xe_R-KPu*Pu_i-KU*U_i)*exp(-t./tau_D)+Xe_R...
        +KPu*Pu_i*exp(-t./tau_Pu)+KU*U_i*exp(-t./tau_U);

Xe_R=Xe131_R1;Xe_i=Xe131_i;
KPu=KPu131;KU=KU131;
Xe131=(Xe_i-Xe_R-KPu*Pu_i-KU*U_i)*exp(-t./tau_D)+Xe_R...
        +KPu*Pu_i*exp(-t./tau_Pu)+KU*U_i*exp(-t./tau_U);

Xe_R=Xe132_R1;Xe_i=Xe132_i;
KPu=KPu132;KU=KU132;
Xe132=(Xe_i-Xe_R-KPu*Pu_i-KU*U_i)*exp(-t./tau_D)+Xe_R...
        +KPu*Pu_i*exp(-t./tau_Pu)+KU*U_i*exp(-t./tau_U);

Xe_R=Xe134_R1;Xe_i=Xe134_i;
KPu=KPu134;KU=KU134;
Xe134=(Xe_i-Xe_R-KPu*Pu_i-KU*U_i)*exp(-t./tau_D)+Xe_R...
        +KPu*Pu_i*exp(-t./tau_Pu)+KU*U_i*exp(-t./tau_U);

Xe_R=Xe136_R1;Xe_i=Xe136_i;
KPu=KPu136;KU=KU136;
Xe136=(Xe_i-Xe_R-KPu*Pu_i-KU*U_i)*exp(-t./tau_D)+Xe_R...
        +KPu*Pu_i*exp(-t./tau_Pu)+KU*U_i*exp(-t./tau_U);

%% highD

% give solution result
% D= % 4~8e24 kg/Ga
% R= % (0.8~0.9)*Xe_pd/Xe_slab*D
% for 128Xe/130Xe, larger D should refer to smaller R
D=10e24; % 5.128
R=0.85*Xe_obs.Xe(2)/Xe_slab*D;
%
% U=
U238bse=13e-9; % [ppb]->[g/g]
U_i=U238bse*exp(Xe_input.lamUr*t_pd/yr_s)/238*NA; % [g/g]->[atmos/g]
Pu_i=U_i*Xe_input.R0Pu_Ur; % [atmos/g]


% calculate dependent params
Xe130_R=Xe_slab*R/D;
% pd atm
Xe128_R1=Xe130_R*Xerpa;
Xe131_R1=Xe130_R*Xe131rpa;
Xe132_R1=Xe130_R*Xe132rpa;
Xe134_R1=Xe130_R*Xe134rpa;
Xe136_R1=Xe130_R*Xe136rpa;
% initial atm
Xe128_R2=Xe130_R*Xe_input.Xes_init;
Xe131_R2=Xe130_R*Xe_input.Xe131r130atm_init;
Xe132_R2=Xe130_R*Xe_input.Xe132r130atm_init;
Xe134_R2=Xe130_R*Xe_input.Xe134r130atm_init;
Xe136_R2=Xe130_R*Xe_input.Xe136r130atm_init;

Xe130_i=Xe_input.Xe_init;
Xe128_i=Xe_input.XeAtm_init;
Xe131_i=Xe_input.Xe131_init;
Xe132_i=Xe_input.Xe132_init;
Xe134_i=Xe_input.Xe134_init;
Xe136_i=Xe_input.Xe136_init;


tau_D=M/D; % [Ga]

KPu130=0;
KU130=0;
KPu128=0;
KU128=0;
KPu131=YPu131/(tau_Pu/tau_D-1);
KU131=YUr131/(tau_U/tau_D-1);
KPu132=YPu132/(tau_Pu/tau_D-1);
KU132=YUr132/(tau_U/tau_D-1);
KPu134=YPu134/(tau_Pu/tau_D-1);
KU134=YUr134/(tau_U/tau_D-1);
KPu136=YPu/(tau_Pu/tau_D-1);
KU136=YUr/(tau_U/tau_D-1);



t=0:0.01:4.6;
% tau_D0=tau_D;
% tau_D=tau_D0/1.5;

Xe_R=Xe130_R;Xe_i=Xe130_i;
KPu=KPu130;KU=KU130;
Xe130_highD=(Xe_i-Xe_R-KPu*Pu_i-KU*U_i)*exp(-t./tau_D)+Xe_R...
        +KPu*Pu_i*exp(-t./tau_Pu)+KU*U_i*exp(-t./tau_U);

Xe_R=Xe128_R1;Xe_i=Xe128_i;
KPu=KPu128;KU=KU128;
Xe128_highD=(Xe_i-Xe_R-KPu*Pu_i-KU*U_i)*exp(-t./tau_D)+Xe_R...
        +KPu*Pu_i*exp(-t./tau_Pu)+KU*U_i*exp(-t./tau_U);


Xe_R=Xe132_R1;Xe_i=Xe132_i;
KPu=KPu132;KU=KU132;
Xe132_highD=(Xe_i-Xe_R-KPu*Pu_i-KU*U_i)*exp(-t./tau_D)+Xe_R...
        +KPu*Pu_i*exp(-t./tau_Pu)+KU*U_i*exp(-t./tau_U);


figure;
set(gcf,'color','white','Units','centimeters','Position',[2 5 17.8 12]);

gca=subplot(2,2,1);
plot(t,Xe130,'k','LineWidth',1.5);
hold on;
plot(t,Xe130_highD,'r','LineWidth',1.5);
% legend(['\tau_D=0.8 Gyr',newline,'Xe_R=1.7e6 atoms/g',newline,'U_{BSE}=13 ppb']);
legend('Analytical Solution','Median Numerical Solution');
ylabel('^{130}Xe Concentration (atoms/g)');
set(gca,'FontUnits','points','FontSize',7,'FontWeight','bold',...
    'LineWidth',1.5,'LabelFontSizeMultiplier',1.0);
text(-0.126,1.01,'a','Units','normalized','FontWeight','bold','FontSize',8);
xlim([0 4.6]);

gca=subplot(2,2,2);
plot(t,Xe132,'k','LineWidth',1.5);
hold on;
plot(t,Xe132_highD,'r','LineWidth',1.5);
% legend(['\tau_D=0.8 Gyr',newline,'Xe_R=1.7e6 atoms/g',newline,'U_{BSE}=13 ppb']);
legend('Analytical Solution','Median Numerical Solution');
ylabel('^{130}Xe Concentration (atoms/g)');
set(gca,'FontUnits','points','FontSize',7,'FontWeight','bold',...
    'LineWidth',1.5,'LabelFontSizeMultiplier',1.0);
text(-0.126,1.01,'a','Units','normalized','FontWeight','bold','FontSize',8);
xlim([0 4.6]);

gca=subplot(2,2,3);
plot(t,Xe128./Xe130,'k','LineWidth',1.5);
hold on;
plot(t,Xe128_highD./Xe130_highD,'r','LineWidth',1.5);
% legend(['\tau_D=0.8 Gyr',newline,'Xe_R=1.7e6 atoms/g',newline,'U_{BSE}=13 ppb']);
legend('Analytical Solution','Median Numerical Solution');
ylabel('^{130}Xe Concentration (atoms/g)');
set(gca,'FontUnits','points','FontSize',7,'FontWeight','bold',...
    'LineWidth',1.5,'LabelFontSizeMultiplier',1.0);
text(-0.126,1.01,'a','Units','normalized','FontWeight','bold','FontSize',8);
xlim([0 4.6]);

gca=subplot(2,2,4);
plot(t,Xe132./Xe130,'k','LineWidth',1.5);
hold on;
plot(t,Xe132_highD./Xe130_highD,'r','LineWidth',1.5);
ylabel('^{130}Xe Concentration (atoms/g)');
set(gca,'FontUnits','points','FontSize',7,'FontWeight','bold',...
    'LineWidth',1.5,'LabelFontSizeMultiplier',1.0);
text(-0.126,1.01,'a','Units','normalized','FontWeight','bold','FontSize',8);
xlim([0 4.6]);








