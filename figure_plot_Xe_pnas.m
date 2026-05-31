%% % figure plot for Xe =pnas=  paper

%% figure #1: fit observations
% 2*3
% cc-formation age, cc-surface age, mantle potential temperature,
% Xe concentration, Xe ratio1, Xe ratio2

% load observations
% cc
par_solver_new;
tmax = t_pd/yr_s/1e9;
dt = 0.001;% length of each timestep
t = 0:dt:tmax;
nt=length(t);
data_formationage = readmatrix(...
    'cc/observation/cc_age.xlsx','Sheet','korenaga18a_Tunmix_orig');
% Surface age distribution data from Roberts & Spencer (2015)
data_zircon_surf = readmatrix(...
    'cc/observation/cc_age.xlsx','Sheet','korenaga18a_T_U_Pb');
% use load_FandS_fun function to set data in the same dimension as the time series
[F_Jun_same,S_Jun_same] = load_FandS_fun(t,nt,data_formationage,data_zircon_surf);

% thermal
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

% Xe
default_input_for_degassing;
% observation of Xe: present day mantle
Xe_obs = struct( 'Xe',       [4.3e5	9.2e5], ... % atoms/g
                 'Xe128v130',     [0.475	0.478], ...
                 'Xe128v132',     [0.069	0.071], ...
                 'Xe130v132',     [0.1445	0.1493], ...
                 'Xe131v132',     [0.7608	0.7786], ...
                 'Xe134v132',     [0.4082	0.4302], ...
                 'Xe136v132',     [0.3559	0.3835] ...
                );
% observation of Xe: starting material of mantle
Xe_start_cc = struct('Xe',       [3.2e7 3.2e8], ... % atoms/g
                 'Xe128v130',     Xercc, ...
                 'Xe128v132',     Xercc./Xe132rcc, ...
                 'Xe130v132',     1./Xe132rcc, ...
                 'Xe131v132',     Xe131rcc./Xe132rcc, ...
                 'Xe134v132',     Xe134rcc./Xe132rcc, ...
                 'Xe136v132',     Xe136rcc./Xe132rcc ...
                );

Xe_start_en = struct('Xe',       [3.2e7 3.2e8], ... % atoms/g
                 'Xe128v130',     Xeren, ...
                 'Xe128v132',     Xeren./Xe132ren, ...
                 'Xe130v132',     1./Xe132ren, ...
                 'Xe131v132',     Xe131ren./Xe132ren, ...
                 'Xe134v132',     Xe134ren./Xe132ren, ...
                 'Xe136v132',     Xe136ren./Xe132ren ...
                );


% observation of Xe: starting material of atmosphere
Xe_atm_start = struct(...
                 'Xe128v130',     Xe_input.Xes_init, ...
                 'Xe128v132',     Xe_input.Xes_init./Xe_input.Xe132r130atm_init, ...
                 'Xe130v132',     1./Xe_input.Xe132r130atm_init, ...
                 'Xe131v132',     Xe_input.Xe131r130atm_init./Xe_input.Xe132r130atm_init, ...
                 'Xe134v132',     Xe_input.Xe134r130atm_init./Xe_input.Xe132r130atm_init, ...
                 'Xe136v132',     Xe_input.Xe136r130atm_init./Xe_input.Xe132r130atm_init ...
                );

% observation of Xe: present day atmosphere
Xe_atm_obs = struct(...
                 'Xe128v130',     Xerpa, ...
                 'Xe128v132',     Xerpa./Xe132rpa, ...
                 'Xe130v132',     1./Xe132rpa, ...
                 'Xe131v132',     Xe131rpa./Xe132rpa, ...
                 'Xe134v132',     Xe134rpa./Xe132rpa, ...
                 'Xe136v132',     Xe136rpa./Xe132rpa ...
                );


% load data
cc_output_margin_struct=load('cc/cc_output_for_figure_margin.mat');
cc_output_margin=cc_output_margin_struct.cc_output_margin;

% thermal_output_margin_struct=load('thermal/thermal_output_for_figure_margin.mat');
% thermal_output_margin=thermal_output_margin_struct.thermal_output_margin;

Xe_output_margin_struct=load('Xe/Xe_output_for_figure_margin.mat');
Xe_output_margin=Xe_output_margin_struct.Xe_output_margin;

Xe_output_margin_struct0=load('Xe/Xe_output_for_figure_margin0.mat');
Xe_output_margin_struct1=load('Xe/Xe_output_for_figure_margin1.mat');

Xe_output_ex_struct=load('Xe/Xe_output_for_figure_ex.mat');
Xe_output_ex=Xe_output_ex_struct.Xe_output;


% cc
F_lower = cc_output_margin.F_Xe_lower;
F_upper = cc_output_margin.F_Xe_upper;
F_5 = cc_output_margin.F_Xe_5;
F_25 = cc_output_margin.F_Xe_25;
F_75 = cc_output_margin.F_Xe_75;
F_95 = cc_output_margin.F_Xe_95;

S_lower = cc_output_margin.S_Xe_lower;
S_upper = cc_output_margin.S_Xe_upper;
S_5 = cc_output_margin.S_Xe_5;
S_25 = cc_output_margin.S_Xe_25;
S_75 = cc_output_margin.S_Xe_75;
S_95 = cc_output_margin.S_Xe_95;


% thermal
Tp_lower = Xe_output_margin_struct0.Xe_output_margin.Tp_Xe_lower;
Tp_upper = Xe_output_margin_struct0.Xe_output_margin.Tp_Xe_upper;
Tp_5 = Xe_output_margin.Tp_Xe_5;
Tp_25 = Xe_output_margin.Tp_Xe_25;
Tp_75 = Xe_output_margin.Tp_Xe_75;
Tp_95 = Xe_output_margin.Tp_Xe_95;

%
% Xe 
Xe_lower = Xe_output_margin_struct1.Xe_output_margin.Xe_lower;
Xe_upper = Xe_output_margin_struct1.Xe_output_margin.Xe_upper;
Xe_5 = Xe_output_margin.Xe_5;
Xe_25 = Xe_output_margin.Xe_25;
Xe_75 = Xe_output_margin.Xe_75;
Xe_95 = Xe_output_margin.Xe_95;

Xe128t=Xe_output_ex.Xe128(:,1);
Xe130t=Xe_output_ex.Xe(:,1);
Xe132t=Xe_output_ex.Xe132(:,1);
Xe136t=Xe_output_ex.Xe136(:,1);

Xe128t_lowcc_kg=Xe_output_ex.Xe128(:,2);
Xe130t_lowcc_kg=Xe_output_ex.Xe(:,2);
Xe132t_lowcc_kg=Xe_output_ex.Xe132(:,2);
Xe136t_lowcc_kg=Xe_output_ex.Xe136(:,2);

Xe128t_lowcc_kr=Xe_output_ex.Xe128(:,3);
Xe130t_lowcc_kr=Xe_output_ex.Xe(:,3);
Xe132t_lowcc_kr=Xe_output_ex.Xe132(:,3);
Xe136t_lowcc_kr=Xe_output_ex.Xe136(:,3);

Xe128t_lowcc_Rs=Xe_output_ex.Xe128(:,4);
Xe130t_lowcc_Rs=Xe_output_ex.Xe(:,4);
Xe132t_lowcc_Rs=Xe_output_ex.Xe132(:,4);
Xe136t_lowcc_Rs=Xe_output_ex.Xe136(:,4);

Xe128t_highU=Xe_output_ex.Xe128(:,5);
Xe130t_highU=Xe_output_ex.Xe(:,5);
Xe132t_highU=Xe_output_ex.Xe132(:,5);
Xe136t_highU=Xe_output_ex.Xe136(:,5);


%
figure;
set(gcf,'color','white');

color_unit=[0.6, 0.8, 1.0];
color_purple = [0.5, 0.2, 0.7];


% cc Formation age
dt = 0.001;% length of each timestep
t = 0:dt:tmax;
% F_5 = cc_output_margin.F_cc_5;
% F_95 = cc_output_margin.F_cc_95;
F_5 = cc_output_margin.F_Xe_5;
F_95 = cc_output_margin.F_Xe_95;
% t = t';
% t = t(:)'; 
subplot(2,3,1);                       
x_fill = [t, fliplr(t)];
y_fill0 = [F_lower, fliplr(F_upper)];
y_fill1 = [F_5, fliplr(F_95)];
y_fill2 = [F_25, fliplr(F_75)];

fill(x_fill, y_fill0, color_unit, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on;
fill(x_fill, y_fill1, color_unit, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on;
fill(x_fill, y_fill2, color_unit, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on;
plot(t, F_Jun_same,'r-');
xlabel('Time (Gyr)');
ylabel('Formation age distribution');
grid on;
xlim([0 4.6]);ylim([0 1.0]);


% surface age
subplot(2,3,2);
x_fill = [t, fliplr(t)];
y_fill0 = [S_lower, fliplr(S_upper)];
y_fill1 = [S_5, fliplr(S_95)];
y_fill2 = [S_25, fliplr(S_75)];
fill(x_fill, y_fill0, color_unit, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on;
fill(x_fill, y_fill1, color_unit, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on;
fill(x_fill, y_fill2, color_unit, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on;
plot(t, S_Jun_same,'r-');
xlabel('Time (Gyr)');
ylabel('Surface age distribution');
grid on;
xlim([0 4.6]);ylim([0 1.0]);


% potential temperature
default_input_for_degassing;
t=Simulation_input.time_series;
t=t/yr_s/1e9;% to Ga
t=t';
subplot(2,3,3);
x_fill = [t, fliplr(t)];
y_fill0 = [Tp_lower,(fliplr(Tp_upper))];
y_fill1 = [Tp_5,(fliplr(Tp_95))];
y_fill2 = [Tp_25,(fliplr(Tp_75))];
fill(x_fill, y_fill0,color_purple, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on;
fill(x_fill, y_fill1,color_purple, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on;
fill(x_fill, y_fill2,color_purple, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on;

% plot(t, Tp_50, 'm-', 'LineWidth', 1.5);
scatter( t_Herz,Tp_Herz, 100, ...
    'o', 'MarkerEdgeColor', [1, 0.5, 0], ...
    'MarkerFaceColor', [1, 0.3, 0], 'LineWidth', 1.5);
% % 
% stairs(t, F_50, 'Color', [1 0.4 0.2], 'LineWidth', 2); 
%
xlabel('Time (Gyr)');
ylabel('Mantle Potential Temperature (°C)');
grid on;
xlim([0 4.6]);


% Xe concentration

subplot(2,3,4);
x_fill = [t, fliplr(t)];
% y_fill0 = [Xe_lower,(fliplr(Xe_upper))];
y_fill1 = [Xe_5,(fliplr(Xe_95))];
y_fill2 = [Xe_25,(fliplr(Xe_75))];

% fill(x_fill, y_fill0,color_purple, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on;
fill(x_fill, y_fill1,color_purple, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on;
fill(x_fill, y_fill2,color_purple, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on;


% Xe ratio1: 128Xe/132Xe vs 130Xe/132Xe
time=tmax-t;

x1=Xe128t./Xe132t;
x2=Xe130t./Xe132t;


x1_lowcc=Xe128t_lowcc_Rs./Xe132t_lowcc_Rs;
x2_lowcc=Xe130t_lowcc_Rs./Xe132t_lowcc_Rs;

x1_highU=Xe128t_highU./Xe132t_highU;
x2_highU=Xe130t_highU./Xe132t_highU;


% present day mantle: a range
x1obs=Xe_obs.Xe128v132;x2obs=Xe_obs.Xe130v132;
% starting mantle: a data ponit
x1start=Xe_start_en.Xe128v132;x2start=Xe_start_en.Xe130v132;
% present day atmosphere: a data point
x1atm_obs=Xe_atm_obs.Xe128v132;x2atm_obs=Xe_atm_obs.Xe130v132;
% starting atmosphere: a data point
x1atm_start=Xe_atm_start.Xe128v132;x2atm_start=Xe_atm_start.Xe130v132;

subplot(2,3,5);
% evolutionary curve for successful params
surface([x1(:) x1(:)], [x2(:) x2(:)], ...
    [zeros(size(x1(:))) zeros(size(x2(:)))], ...
        [time(:) time(:)], 'EdgeColor', ...
        'interp', 'FaceColor', 'none', 'LineWidth', 2);
hold on;
% evolutionary curve for low cc params
surface([x1_lowcc(:) x1_lowcc(:)], [x2_lowcc(:) x2_lowcc(:)], ...
    [zeros(size(x1_lowcc(:))) zeros(size(x2_lowcc(:)))], ...
        [time(:) time(:)], 'EdgeColor', ...
        'interp', 'FaceColor', 'none', 'LineWidth', 2);
hold on;
% evolutionary curve for high Ubse params
surface([x1_highU(:) x1_highU(:)], [x2_highU(:) x2_highU(:)], ...
    [zeros(size(x1_highU(:))) zeros(size(x2_highU(:)))], ...
        [time(:) time(:)], 'EdgeColor', ...
        'interp', 'FaceColor', 'none', 'LineWidth', 2);
hold on;

% present day mantle
fill([x1obs(1) x1obs(2) x1obs(2) x1obs(1)], [x2obs(1) x2obs(1) x2obs(2) x2obs(2)], [0.5 0.5 0.5], ...
     'FaceAlpha', 0.5, 'EdgeColor', 'none');hold on;
% starting mantle
scatter(x1start, x2start, 150, 'o', 'MarkerFaceColor',[0.2 0.6 0.8], ...
        'MarkerEdgeColor','k', 'LineWidth',1.2);hold on;
% present day atmosphere
scatter(x1atm_obs, x2atm_obs, 150, 'o', 'MarkerFaceColor',[0.2 0.6 0.8], ...
        'MarkerEdgeColor','k', 'LineWidth',1.2);hold on;
% starting mantle
scatter(x1atm_start, x2atm_start, 150, 'o', 'MarkerFaceColor',[0.2 0.6 0.8], ...
        'MarkerEdgeColor','k', 'LineWidth',1.2);hold on;


xlim([0.060 0.085]);ylim([0.12 0.17]);
colormap(jet);
colorbar;grid on;
xlabel('^{128}Xe/^{132}Xe');
ylabel('^{130}Xe/^{132}Xe');

%
% Xe ratio2: 128Xe/132Xe vs 136Xe/132Xe
time=tmax-t;

x1=Xe128t./Xe132t;
x2=Xe136t./Xe132t;


x1_lowcc=Xe128t_lowcc_Rs./Xe132t_lowcc_Rs;
x2_lowcc=Xe136t_lowcc_Rs./Xe132t_lowcc_Rs;

x1_highU=Xe128t_highU./Xe132t_highU;
x2_highU=Xe136t_highU./Xe132t_highU;


% present day mantle: a range
x1obs=Xe_obs.Xe128v132;x2obs=Xe_obs.Xe136v132;
% starting mantle: a data ponit
x1start=Xe_start_en.Xe128v132;x2start=Xe_start_en.Xe136v132;
% present day atmosphere: a data point
x1atm_obs=Xe_atm_obs.Xe128v132;x2atm_obs=Xe_atm_obs.Xe136v132;
% starting atmosphere: a data point
x1atm_start=Xe_atm_start.Xe128v132;x2atm_start=Xe_atm_start.Xe136v132;

subplot(2,3,6);
% evolutionary curve for successful params
surface([x1(:) x1(:)], [x2(:) x2(:)], ...
    [zeros(size(x1(:))) zeros(size(x2(:)))], ...
        [time(:) time(:)], 'EdgeColor', ...
        'interp', 'FaceColor', 'none', 'LineWidth', 2);
hold on;
% evolutionary curve for low cc params
surface([x1_lowcc(:) x1_lowcc(:)], [x2_lowcc(:) x2_lowcc(:)], ...
    [zeros(size(x1_lowcc(:))) zeros(size(x2_lowcc(:)))], ...
        [time(:) time(:)], 'EdgeColor', ...
        'interp', 'FaceColor', 'none', 'LineWidth', 2);
hold on;
% evolutionary curve for high Ubse params
surface([x1_highU(:) x1_highU(:)], [x2_highU(:) x2_highU(:)], ...
    [zeros(size(x1_highU(:))) zeros(size(x2_highU(:)))], ...
        [time(:) time(:)], 'EdgeColor', ...
        'interp', 'FaceColor', 'none', 'LineWidth', 2);
hold on;

% present day mantle
fill([x1obs(1) x1obs(2) x1obs(2) x1obs(1)], [x2obs(1) x2obs(1) x2obs(2) x2obs(2)], [0.5 0.5 0.5], ...
     'FaceAlpha', 0.5, 'EdgeColor', 'none');hold on;
% starting mantle
scatter(x1start, x2start, 150, 'o', 'MarkerFaceColor',[0.2 0.6 0.8], ...
        'MarkerEdgeColor','k', 'LineWidth',1.2);hold on;
% present day atmosphere
scatter(x1atm_obs, x2atm_obs, 150, 'o', 'MarkerFaceColor',[0.2 0.6 0.8], ...
        'MarkerEdgeColor','k', 'LineWidth',1.2);hold on;
% starting mantle
scatter(x1atm_start, x2atm_start, 150, 'o', 'MarkerFaceColor',[0.2 0.6 0.8], ...
        'MarkerEdgeColor','k', 'LineWidth',1.2);hold on;


xlim([0.060 0.085]);ylim([0.25 0.65]);
colormap(jet);
colorbar;grid on;
xlabel('^{128}Xe/^{132}Xe');
ylabel('^{136}Xe/^{132}Xe');


%% figure #2: importance & distribution of params
% clc;clear;

% load data
% Xe_rf_struct=load('Xe/Xe_rf_residuals_Xei3p2e7_alot.mat');
Xe_rf_importance=Xe_rf_struct.RF_importance;
Xe_rf_target=Xe_rf_struct.resid_col_names;
% Xe_rf_params=Xe_rf_struct.param_names;

Xe_rf_params= {'t_{s}','\kappa_g','\kappa_r','R_s','R_p',...
    'U_{bse}','Qc', 'Ti', '\eta_{ref}',...
    'Fr',...
    'Fd_{M}','Fd_{P}'};

% load data
success_F_struct=load('cc/out_cc_grwoth_success_lhs_MC_Krw.mat');
success_F=success_F_struct.success_F;
[size_F,nparam_cc]=size(success_F);

success_T_struct=load('thermal/out_thermal_success_lhs_test_fixK_U_MC_lowerT.mat');
success_T=success_T_struct.success_T;
[size_T,nparam_T]=size(success_T);

success_Xe_struct=load('Xe/out_Xe_success_lhs_new_3e7_test_fixK_U_MC_negativFD.mat');
success_Xe=success_Xe_struct.success_Xe;
[size_Xe,nparam_Xe]=size(success_Xe);


%%
figure;
set(gcf,'color','white');

% Feature importance of ^{130}Xe l=1
subplot(2,3,1);
l=1;
bar(Xe_rf_importance(l, :)); 
xticks(1:length(Xe_rf_params));
set(gca, 'XTickLabel', Xe_rf_params);
ylabel('Feature Importance');
title(['Influence on ', Xe_rf_target{l}], 'FontSize', 16,'FontWeight', 'bold');
grid on;

% Feature importance of ^{128}Xe/^{130}Xe l=2
subplot(2,3,2);
l=2;
bar(Xe_rf_importance(l, :)); 
xticks(1:length(Xe_rf_params));
set(gca, 'XTickLabel', Xe_rf_params);
ylabel('Feature Importance');
title(['Influence on ', Xe_rf_target{l}], 'FontSize', 16,'FontWeight', 'bold');
grid on;

% Feature importance of '^{130}Xe/^{132}Xe' l=3
subplot(2,3,3);
l=3;
bar(Xe_rf_importance(l, :)); 
xticks(1:length(Xe_rf_params));
set(gca, 'XTickLabel', Xe_rf_params);
ylabel('Feature Importance');
title(['Influence on ', Xe_rf_target{l}], 'FontSize', 16,'FontWeight', 'bold');
grid on;

% histogram of Ubse
k=6;
y2=success_T(:,k);
y3=success_Xe(:,k);

subplot(2,3,4);
histogram(y2, 30, 'Normalization','probability',...            
              'EdgeColor','b', ...      % boundary
              'FaceColor','none', ...   % fill
              'LineWidth',1.2);         % boundary
hold on;% 'BinWidth',0.1, ...
histogram(y3, 30, 'Normalization','probability',...
              'EdgeColor','r', ...      % boundary
              'FaceColor','none', ...   % fill
              'LineWidth',1.2);         % boundary


% histogram of kappa_rcc
k=3;
y1=success_F(:,k);
y2=success_T(:,k);
y3=success_Xe(:,k);

subplot(2,3,5);
histogram(y1, 30, 'Normalization','probability',...            
              'EdgeColor','k', ...      % boundary
              'FaceColor','none', ...   % fill
              'LineWidth',1.2);         % boundary
hold on;% 'BinWidth',0.1, ...
histogram(y2, 30, 'Normalization','probability',...
              'EdgeColor','b', ...      % boundary
              'FaceColor','none', ...   % fill
              'LineWidth',1.2);         % boundary
histogram(y3, 30, 'Normalization','probability',...
              'EdgeColor','r', ...      % boundary
              'FaceColor','none', ...   % fill
              'LineWidth',1.2);         % boundary



% histogram of Rs
k=4;
y1=success_F(:,k);
y2=success_T(:,k);
y3=success_Xe(:,k);

subplot(2,3,6);
histogram(y1, 30, 'Normalization','probability',...            
              'EdgeColor','k', ...      % boundary
              'FaceColor','none', ...   % fill
              'LineWidth',1.2);         % boundary
hold on;% 'BinWidth',0.1, ...
histogram(y2, 30, 'Normalization','probability',...
              'EdgeColor','b', ...      % boundary
              'FaceColor','none', ...   % fill
              'LineWidth',1.2);         % boundary
histogram(y3, 30, 'Normalization','probability',...
              'EdgeColor','r', ...      % boundary
              'FaceColor','none', ...   % fill
              'LineWidth',1.2);         % boundary


%% figure #3: cc curves after vs before; degassing analysis
% clc;clear;
par_solver_new;
tmax = t_pd/yr_s/1e9;
dt = 0.001;% length of each timestep
t = 0:dt:tmax;
nt=length(t);

% load data
cc_output_margin_struct=load('cc/cc_output_for_figure_margin.mat');
cc_output_margin=cc_output_margin_struct.cc_output_margin;

% thermal_output_margin_struct=load('thermal/thermal_output_for_figure_margin.mat');
% thermal_output_margin=thermal_output_margin_struct.thermal_output_margin;

Xe_output_margin_struct=load('Xe/Xe_output_for_figure_margin.mat');
Xe_output_margin=Xe_output_margin_struct.Xe_output_margin;

Xe_output_margin_struct0=load('Xe/Xe_output_for_figure_margin0.mat');
Xe_output_margin_struct1=load('Xe/Xe_output_for_figure_margin1.mat');

Xe_output_ex_struct=load('Xe/Xe_output_for_figure_ex.mat');
Xe_output_ex=Xe_output_ex_struct.Xe_output;

%
NGcc_cc_lower = cc_output_margin.NGcc_lower;
NGcc_cc_upper = cc_output_margin.NGcc_upper;
NGcc_cc_5 = cc_output_margin.NGcc_cc_5;
NGcc_cc_25 = cc_output_margin.NGcc_cc_25;
NGcc_cc_75 = cc_output_margin.NGcc_cc_75;
NGcc_cc_95 = cc_output_margin.NGcc_cc_95;

Gcc_cc_lower = cc_output_margin.Gcc_lower;
Gcc_cc_upper = cc_output_margin.Gcc_upper;
Gcc_cc_5 = cc_output_margin.Gcc_cc_5;
Gcc_cc_25 = cc_output_margin.Gcc_cc_25;
Gcc_cc_75 = cc_output_margin.Gcc_cc_75;
Gcc_cc_95 = cc_output_margin.Gcc_cc_95;

Rcc_cc_lower = cc_output_margin.Rcc_lower;
Rcc_cc_upper = cc_output_margin.Rcc_upper;
Rcc_cc_5 = cc_output_margin.Rcc_cc_5;
Rcc_cc_25 = cc_output_margin.Rcc_cc_25;
Rcc_cc_75 = cc_output_margin.Rcc_cc_75;
Rcc_cc_95 = cc_output_margin.Rcc_cc_95;

NGcc_Xe_lower = cc_output_margin.NGcc_Xe_lower;
NGcc_Xe_upper = cc_output_margin.NGcc_Xe_upper;
NGcc_Xe_5 = cc_output_margin.NGcc0_Xe_5;
NGcc_Xe_25 = cc_output_margin.NGcc0_Xe_25;
NGcc_Xe_75 = cc_output_margin.NGcc0_Xe_75;
NGcc_Xe_95 = cc_output_margin.NGcc0_Xe_95;

Gcc_Xe_lower = cc_output_margin.Gcc_Xe_lower;
Gcc_Xe_upper = cc_output_margin.Gcc_Xe_upper;
Gcc_Xe_5 = cc_output_margin.Gcc0_Xe_5;
Gcc_Xe_25 = cc_output_margin.Gcc0_Xe_25;
Gcc_Xe_75 = cc_output_margin.Gcc0_Xe_75;
Gcc_Xe_95 = cc_output_margin.Gcc0_Xe_95;

Rcc_Xe_lower = cc_output_margin.Rcc_Xe_lower;
Rcc_Xe_upper = cc_output_margin.Rcc_Xe_upper;
Rcc_Xe_5 = cc_output_margin.Rcc0_Xe_5;
Rcc_Xe_25 = cc_output_margin.Rcc0_Xe_25;
Rcc_Xe_75 = cc_output_margin.Rcc0_Xe_75;
Rcc_Xe_95 = cc_output_margin.Rcc0_Xe_95;

%
Dcc_lower = Xe_output_margin.D0_cc_lower;
Dcc_upper = Xe_output_margin.D0_cc_upper;
Dcc_5 = Xe_output_margin.D0_cc_5;
Dcc_25 = Xe_output_margin.D0_cc_25;
Dcc_50 = Xe_output_margin.D0_cc_50;
Dcc_75 = Xe_output_margin.D0_cc_75;
Dcc_95 = Xe_output_margin.D0_cc_95;

Dmor_lower = Xe_output_margin.D0_mor_lower;
Dmor_upper = Xe_output_margin.D0_mor_upper;
Dmor_5 = Xe_output_margin.D0_mor_5;
Dmor_50 = Xe_output_margin.D0_mor_50;
Dmor_25 = Xe_output_margin.D0_mor_25;
Dmor_75 = Xe_output_margin.D0_mor_75;
Dmor_95 = Xe_output_margin.D0_mor_95;

Dplume_lower = Xe_output_margin.D0_p_lower;
Dplume_upper = Xe_output_margin.D0_p_upper;
Dplume_5 = Xe_output_margin.D0_p_5;
Dplume_50 = Xe_output_margin.D0_p_50;
Dplume_25 = Xe_output_margin.D0_p_25;
Dplume_75 = Xe_output_margin.D0_p_75;
Dplume_95 = Xe_output_margin.D0_p_95;




%%
figure;
set(gcf,'color','white');
color_unit=[0.6, 0.8, 1.0];
color_purple = [0.5, 0.2, 0.7];
%% % Gcc
%subplot(2,2,1);
x_fill = [t, fliplr(t)];

y_fill0 = [Gcc_Xe_lower,(fliplr(Gcc_Xe_upper))];
y_fill1 = [Gcc_Xe_5,(fliplr(Gcc_Xe_95))];
y_fill2 = [Gcc_Xe_25,(fliplr(Gcc_Xe_75))];

z_fill0 = [Gcc_cc_lower,(fliplr(Gcc_cc_upper))];
z_fill1 = [Gcc_cc_5,(fliplr(Gcc_cc_95))];
z_fill2 = [Gcc_cc_25,(fliplr(Gcc_cc_75))];

fill(x_fill, z_fill0,color_unit, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on;
% fill(x_fill, z_fill1,color_unit, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on;
% fill(x_fill, z_fill2,color_unit, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on;

 fill(x_fill, y_fill0,color_purple, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on;
% fill(x_fill, y_fill1,color_purple, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on;
% fill(x_fill, y_fill2,color_purple, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on;
xlabel('Time (Gyr)'); ylabel('Crustal generation rate (kg/Ga)'); grid on;



%% % Rcc
%subplot(2,2,2);
x_fill = [t, fliplr(t)];

y_fill0 = [Rcc_Xe_lower,(fliplr(Rcc_Xe_upper))];
y_fill1 = [Rcc_Xe_5,(fliplr(Rcc_Xe_95))];
y_fill2 = [Rcc_Xe_25,(fliplr(Rcc_Xe_75))];

z_fill0 = [Rcc_cc_lower,(fliplr(Rcc_cc_upper))];
z_fill1 = [Rcc_cc_5,(fliplr(Rcc_cc_95))];
z_fill2 = [Rcc_cc_25,(fliplr(Rcc_cc_75))];

 fill(x_fill, z_fill0,color_unit, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on;
% fill(x_fill, z_fill1,color_unit, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on;
% fill(x_fill, z_fill2,color_unit, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on;

 fill(x_fill, y_fill0,color_purple, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on;
% fill(x_fill, y_fill1,color_purple, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on;
% fill(x_fill, y_fill2,color_purple, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on;
xlabel('Time (Gyr)'); ylabel('Crustal recycling rate (kg/Ga)'); grid on;


%% % NGcc
%subplot(2,2,3);
x_fill = [t, fliplr(t)];

y_fill0 = [NGcc_Xe_lower,(fliplr(NGcc_Xe_upper))];
y_fill1 = [NGcc_Xe_5,(fliplr(NGcc_Xe_95))];
y_fill2 = [NGcc_Xe_25,(fliplr(NGcc_Xe_75))];

z_fill0 = [NGcc_cc_lower,(fliplr(NGcc_cc_upper))];
z_fill1 = [NGcc_cc_5,(fliplr(NGcc_cc_95))];
z_fill2 = [NGcc_cc_25,(fliplr(NGcc_cc_75))];

 fill(x_fill, z_fill0,color_unit, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on;
% fill(x_fill, z_fill1,color_unit, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on;
% fill(x_fill, z_fill2,color_unit, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on;

 fill(x_fill, y_fill0,color_purple, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on;
% fill(x_fill, y_fill1,color_purple, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on;
% fill(x_fill, y_fill2,color_purple, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on;
xlabel('Time (Gyr)'); ylabel('Net crustal growth (kg)'); grid on;


%% % Degassing analysis
default_input_for_degassing;
t=Simulation_input.time_series;
t=t/yr_s/1e9;% to Ga
t=t';

% subplot(2,2,4);
plot(t,Dcc_50,'r');hold on;
plot(t,Dmor_50,'b');hold on;
plot(t,Dplume_50,'c');hold on;


%% % figure #4: why highcc and lowU
%% highcc 
% clc;clear all;
contour_struct=load('Xe/fin_6076_KgFr_Xe_3e7.mat');
results=contour_struct.results;
col_names=contour_struct.col_names;
Xq=contour_struct.FrXe;
Yq=contour_struct.Yq;

%
subplot(2,2,1);
Erxec_idx = find(strcmp(col_names, 'Er_xec'));
Erxe128v130_idx = find(strcmp(col_names, 'Er_xe128v130'));
Erxe128v132_idx = find(strcmp(col_names, 'Er_xe128v132'));
Erxe130v132_idx = find(strcmp(col_names, 'Er_xe130v132'));
Erxe131v132_idx = find(strcmp(col_names, 'Er_xe131v132'));
Erxe134v132_idx = find(strcmp(col_names, 'Er_xe134v132'));
Erxe136v132_idx = find(strcmp(col_names, 'Er_xe136v132'));

resid_col_names = {'^{130}Xe','^{128}Xe/^{130}Xe','^{130}Xe/^{132}Xe','^{131}Xe/^{132}Xe','^{134}Xe/^{132}Xe','^{136}Xe/^{132}Xe'};
resid_col_ids = [Erxec_idx,Erxe128v130_idx, Erxe130v132_idx, Erxe131v132_idx, Erxe134v132_idx, Erxe136v132_idx];

%
% Xq=FrXe;
subset1 = results;

Erxec_e = subset1(:, Erxec_idx);
Erxe128v130_e = subset1(:, Erxe128v130_idx);
Erxe130v132_e = subset1(:, Erxe130v132_idx);
Erxe131v132_e = subset1(:, Erxe131v132_idx);
Erxe134v132_e = subset1(:, Erxe134v132_idx);
Erxe136v132_e = subset1(:, Erxe136v132_idx);

 

Erxec_e_grid = reshape(Erxec_e,[11,11]);
Erxe128v130_e_grid = reshape(Erxe128v130_e,[11,11]);
Erxe130v132_e_grid = reshape(Erxe130v132_e,[11,11]);
Erxe131v132_e_grid = reshape(Erxe131v132_e,[11,11]);
Erxe134v132_e_grid = reshape(Erxe134v132_e,[11,11]);
Erxe136v132_e_grid = reshape(Erxe136v132_e,[11,11]);

 


% overlap

color_unit1=[0, 0.4, 0.8];%blue
color_unit2=[0.8, 0.4, 0];%orange
color_unit3=[0.6 0.1 0.1];%red
[ax1,img1]=plot_regionn(Xq, Yq, Erxec_e_grid ,max(-1,min(min(Erxec_e_grid))),min(1,max(max(Erxec_e_grid))),color_unit1, '--'); % 4.3e5,9.2e5
hold on;
[ax2,img2]=plot_regionn(Xq, Yq, Erxe128v130_e_grid, max(-1,min(min(Erxe128v130_e_grid))),min(1,max(max(Erxe128v130_e_grid))),color_unit3, ':');  % 0.475,0.478
hold on;
[ax3,img3]=plot_regionn(Xq, Yq, Erxe130v132_e_grid, max(-1,min(min(Erxe130v132_e_grid))),min(1,max(max(Erxe130v132_e_grid))),color_unit2, '-.'); % 0.1445,0.1493
hold on;
[ax4,img4]=plot_regionn(Xq, Yq, Erxe131v132_e_grid, max(-1,min(min(Erxe131v132_e_grid))),min(1,max(max(Erxe131v132_e_grid))),color_unit2, '--'); % 0.7608,0.7786
hold on;
[ax5,img5]=plot_regionn(Xq, Yq, Erxe134v132_e_grid, max(-1,min(min(Erxe134v132_e_grid))),min(1,max(max(Erxe134v132_e_grid))),color_unit2, ':');  % 0.4082,0.4302
hold on;
[ax6,img6]=plot_regionn(Xq, Yq, Erxe136v132_e_grid,max(-1,min(min(Erxe136v132_e_grid))),min(1,max(max(Erxe136v132_e_grid))), color_unit2, '-.'); % 0.3559,0.3835

% figure;
% set(gcf,'color','w');
ax = axes;
alpha=0.2;
minX=min(min(Xq));maxX=max(max(Xq));
minY=min(min(Yq));maxY=max(max(Yq));

image(ax, img1, 'AlphaData', alpha*2,'XData', [minX,maxX], 'YData', [maxY,minY]); % 设置透明度
hold on;
image(ax, img2, 'AlphaData', alpha,'XData', [minX,maxX], 'YData', [maxY,minY]);
axis(ax, 'tight');
hold on;
image(ax, img3, 'AlphaData', alpha,'XData', [minX,maxX], 'YData', [maxY,minY]);
axis(ax, 'tight');
hold on;
image(ax, img4, 'AlphaData', alpha,'XData', [minX,maxX], 'YData', [maxY,minY]);
axis(ax, 'tight');
hold on;
image(ax, img5, 'AlphaData', alpha,'XData', [minX,maxX], 'YData', [maxY,minY]);
axis(ax, 'tight');
hold on;
image(ax, img6, 'AlphaData', alpha,'XData', [minX,maxX], 'YData', [maxY,minY]);

contour(Xq, Yq, Erxec_e_grid, [max(-1,min(min(Erxec_e_grid))),max(-1,min(min(Erxec_e_grid)))], '--', 'Color', color_unit1, 'LineWidth', 1.5);
contour(Xq, Yq, Erxec_e_grid, [min(1,max(max(Erxec_e_grid))),min(1,max(max(Erxec_e_grid)))], '--', 'Color', color_unit1, 'LineWidth', 1.5);
contour(Xq, Yq, Erxe128v130_e_grid, [max(-1,min(min(Erxe128v130_e_grid))),max(-1,min(min(Erxe128v130_e_grid)))], '-.', 'Color', color_unit3, 'LineWidth', 1.5);
contour(Xq, Yq, Erxe128v130_e_grid, [min(1,max(max(Erxe128v130_e_grid))),min(1,max(max(Erxe128v130_e_grid)))], '-.', 'Color', color_unit3, 'LineWidth', 1.5);


%hold on;
%image(ax, img_overlap, 'AlphaData', alpha,'XData', [0,1], 'YData', [1,0]);
set(gca,'ydir','normal');
axis(ax, 'tight');
xlabel('FrXe', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('\kappa_{gcc}', 'FontSize', 12, 'FontWeight', 'bold')
title('test',...
    'FontSize', 14, 'FontWeight', 'bold')%Xeini=3.2e7atoms/gram, 


%% low U 
% 
contour_struct=load('Xe/fin_6076_U_Fr_3e7.mat');
results=contour_struct.results;
col_names=contour_struct.col_names;
Xq=contour_struct.FrXe;
Yq=contour_struct.Yq;

%
subplot(2,2,1);
Erxec_idx = find(strcmp(col_names, 'Er_xec'));
Erxe128v130_idx = find(strcmp(col_names, 'Er_xe128v130'));
Erxe128v132_idx = find(strcmp(col_names, 'Er_xe128v132'));
Erxe130v132_idx = find(strcmp(col_names, 'Er_xe130v132'));
Erxe131v132_idx = find(strcmp(col_names, 'Er_xe131v132'));
Erxe134v132_idx = find(strcmp(col_names, 'Er_xe134v132'));
Erxe136v132_idx = find(strcmp(col_names, 'Er_xe136v132'));

resid_col_names = {'^{130}Xe','^{128}Xe/^{130}Xe','^{130}Xe/^{132}Xe','^{131}Xe/^{132}Xe','^{134}Xe/^{132}Xe','^{136}Xe/^{132}Xe'};
resid_col_ids = [Erxec_idx,Erxe128v130_idx, Erxe130v132_idx, Erxe131v132_idx, Erxe134v132_idx, Erxe136v132_idx];

%
% Xq=FrXe;
subset1 = results;

Erxec_e = subset1(:, Erxec_idx);
Erxe128v130_e = subset1(:, Erxe128v130_idx);
Erxe130v132_e = subset1(:, Erxe130v132_idx);
Erxe131v132_e = subset1(:, Erxe131v132_idx);
Erxe134v132_e = subset1(:, Erxe134v132_idx);
Erxe136v132_e = subset1(:, Erxe136v132_idx);

 

Erxec_e_grid = reshape(Erxec_e,[11,11]);
Erxe128v130_e_grid = reshape(Erxe128v130_e,[11,11]);
Erxe130v132_e_grid = reshape(Erxe130v132_e,[11,11]);
Erxe131v132_e_grid = reshape(Erxe131v132_e,[11,11]);
Erxe134v132_e_grid = reshape(Erxe134v132_e,[11,11]);
Erxe136v132_e_grid = reshape(Erxe136v132_e,[11,11]);

 


%% overlap

color_unit1=[0, 0.4, 0.8];%blue
color_unit2=[0.8, 0.4, 0];%orange
color_unit3=[0.6 0.1 0.1];%red
[ax1,img1]=plot_regionn(Xq, Yq, Erxec_e_grid ,max(-1,min(min(Erxec_e_grid))),min(1,max(max(Erxec_e_grid))),color_unit1, '--'); % 4.3e5,9.2e5
hold on;
[ax2,img2]=plot_regionn(Xq, Yq, Erxe128v130_e_grid, max(-1,min(min(Erxe128v130_e_grid))),min(1,max(max(Erxe128v130_e_grid))),color_unit3, ':');  % 0.475,0.478
hold on;
[ax3,img3]=plot_regionn(Xq, Yq, Erxe130v132_e_grid, max(-1,min(min(Erxe130v132_e_grid))),min(1,max(max(Erxe130v132_e_grid))),color_unit2, '-.'); % 0.1445,0.1493
hold on;
[ax4,img4]=plot_regionn(Xq, Yq, Erxe131v132_e_grid, max(-1,min(min(Erxe131v132_e_grid))),min(1,max(max(Erxe131v132_e_grid))),color_unit2, '--'); % 0.7608,0.7786
hold on;
[ax5,img5]=plot_regionn(Xq, Yq, Erxe134v132_e_grid, max(-1,min(min(Erxe134v132_e_grid))),min(1,max(max(Erxe134v132_e_grid))),color_unit2, ':');  % 0.4082,0.4302
hold on;
[ax6,img6]=plot_regionn(Xq, Yq, Erxe136v132_e_grid,max(-1,min(min(Erxe136v132_e_grid))),min(1,max(max(Erxe136v132_e_grid))), color_unit2, '-.'); % 0.3559,0.3835

% figure;
% set(gcf,'color','w');
ax = axes;
alpha=0.2;
minX=min(min(Xq));maxX=max(max(Xq));
minY=min(min(Yq));maxY=max(max(Yq));

image(ax, img1, 'AlphaData', alpha*2,'XData', [minX,maxX], 'YData', [maxY,minY]); % 设置透明度
hold on;
image(ax, img2, 'AlphaData', alpha,'XData', [minX,maxX], 'YData', [maxY,minY]);
axis(ax, 'tight');
hold on;
image(ax, img3, 'AlphaData', alpha,'XData', [minX,maxX], 'YData', [maxY,minY]);
axis(ax, 'tight');
hold on;
image(ax, img4, 'AlphaData', alpha,'XData', [minX,maxX], 'YData', [maxY,minY]);
axis(ax, 'tight');
hold on;
image(ax, img5, 'AlphaData', alpha,'XData', [minX,maxX], 'YData', [maxY,minY]);
axis(ax, 'tight');
hold on;
image(ax, img6, 'AlphaData', alpha,'XData', [minX,maxX], 'YData', [maxY,minY]);

contour(Xq, Yq, Erxec_e_grid, [max(-1,min(min(Erxec_e_grid))),max(-1,min(min(Erxec_e_grid)))], '--', 'Color', color_unit1, 'LineWidth', 1.5);
contour(Xq, Yq, Erxec_e_grid, [min(1,max(max(Erxec_e_grid))),min(1,max(max(Erxec_e_grid)))], '--', 'Color', color_unit1, 'LineWidth', 1.5);
contour(Xq, Yq, Erxe128v130_e_grid, [max(-1,min(min(Erxe128v130_e_grid))),max(-1,min(min(Erxe128v130_e_grid)))], '-.', 'Color', color_unit3, 'LineWidth', 1.5);
contour(Xq, Yq, Erxe128v130_e_grid, [min(1,max(max(Erxe128v130_e_grid))),min(1,max(max(Erxe128v130_e_grid)))], '-.', 'Color', color_unit3, 'LineWidth', 1.5);


%hold on;
%image(ax, img_overlap, 'AlphaData', alpha,'XData', [0,1], 'YData', [1,0]);
set(gca,'ydir','normal');
axis(ax, 'tight');
xlabel('FrXe', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('U_{BSE}', 'FontSize', 12, 'FontWeight', 'bold')
title('test',...
    'FontSize', 14, 'FontWeight', 'bold')%Xeini=3.2e7atoms/gram, 






