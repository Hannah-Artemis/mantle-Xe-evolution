%% % figure plot for Xe =pnas=  paper

%% figure S1 &  #1: fit observations
% 2*3

% load observations
% cc

% load extra colormap
load('glasgow.mat');


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

data_Tp2 = readmatrix(...
    'thermal/observation/Tp.xlsx','Sheet','Ganne_Feng');
[Tp_anchorHerz1,Tp_anchorHerz2,Tp_anchorHerz3,Tp_anchorHerz4,...
    t_anchorHerz1,t_anchorHerz2,t_anchorHerz3,t_anchorHerz4,...
    t_GF,Tp_GF] = load_Tp_fun(data_Tp2,tmax); %
t_GF_fin=t_GF(t_GF<4 & data_Tp2(:,5)>10 & data_Tp2(:,6)~=0);
Tp_GF_fin=Tp_GF(t_GF<4 & data_Tp2(:,5)>10 & data_Tp2(:,6)~=0);

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

cc_output_margin_struct=load('Xe/cc_output_for_figure_margin_fin0_lowf.mat');
cc_output_margin=cc_output_margin_struct.cc_output_margin;



Xe_output_margin_struct=load('Xe/Xe_output_for_figure_margin_new_lowf.mat');
Xe_output_margin=Xe_output_margin_struct.Xe_output_margin;


Xe_output_ex_struct=load('Xe/Xe_output_for_figure_ex_new_lowf.mat');
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
Tp_lower = Xe_output_margin_struct.Xe_output_margin.Tp_Xe_lower;
Tp_upper = Xe_output_margin_struct.Xe_output_margin.Tp_Xe_upper;
Tp_5 = Xe_output_margin.Tp_Xe_5;
Tp_25 = Xe_output_margin.Tp_Xe_25;
Tp_75 = Xe_output_margin.Tp_Xe_75;
Tp_95 = Xe_output_margin.Tp_Xe_95;

%
% Xe 
% calc_analytical_solution_ex;
Xe_lower = Xe_output_margin_struct.Xe_output_margin.Xe_lower;
Xe_upper = Xe_output_margin_struct.Xe_output_margin.Xe_upper;
Xe_5 = Xe_output_margin.Xe_5;
Xe_25 = Xe_output_margin.Xe_25;
Xe_75 = Xe_output_margin.Xe_75;
Xe_95 = Xe_output_margin.Xe_95;


Xe128r132_5 = Xe_output_margin.Xe128r132_5;
Xe128r132_25 = Xe_output_margin.Xe128r132_25;
Xe128r132_50 = Xe_output_margin.Xe128r132_50;
Xe128r132_75 = Xe_output_margin.Xe128r132_75;
Xe128r132_95 = Xe_output_margin.Xe128r132_95;

Xe130r132_5 = Xe_output_margin.Xer132_5;
Xe130r132_25 = Xe_output_margin.Xer132_25;
Xe130r132_50 = Xe_output_margin.Xer132_50;
Xe130r132_75 = Xe_output_margin.Xer132_75;
Xe130r132_95 = Xe_output_margin.Xer132_95;

Xe136r132_5 = Xe_output_margin.Xe136r132_5;
Xe136r132_25 = Xe_output_margin.Xe136r132_25;
Xe136r132_50 = Xe_output_margin.Xe136r132_50;
Xe136r132_75 = Xe_output_margin.Xe136r132_75;
Xe136r132_95 = Xe_output_margin.Xe136r132_95;


Xe128t=Xe_output_ex.Xe128(:,1);
Xe130t=Xe_output_ex.Xe(:,1);
Xe132t=Xe_output_ex.Xe132(:,1);
Xe136t=Xe_output_ex.Xe136(:,1);

Xe128r132_atm=Xe_output_ex.Xe128r132_atm;
Xe130r132_atm=Xe_output_ex.Xe130r132_atm;
Xe136r132_atm=Xe_output_ex.Xe136r132_atm;

%% plot sup figure: analytical solution vs median numerical solution

figure;
set(gcf,'color','white','Units','centimeters','Position',[2 5 17.8 12]);

tt=0:0.01:4.6;
gca=subplot(2,2,1);
plot(tt,Xe130,'k','LineWidth',1.5);
hold on;
plot(t_n,Xe_50,'r','LineWidth',1.5);
% legend(['\tau_D=0.8 Gyr',newline,'Xe_R=1.7e6 atoms/g',newline,'U_{BSE}=13 ppb']);
legend('Analytical Solution','Median Numerical Solution');
ylabel('^{130}Xe Concentration (atoms/g)');
set(gca,'FontUnits','points','FontSize',7,'FontWeight','bold',...
    'LineWidth',1.5,'LabelFontSizeMultiplier',1.0);
text(-0.126,1.01,'a','Units','normalized','FontWeight','bold','FontSize',8);
xlim([0 4.6]);

gca=subplot(2,2,3);
plot(tt,Xe128./Xe130,'k','LineWidth',1.5);
hold on;
% plot(t_n,Xe128_50./Xe_50,'k','LineWidth',1.5);
hold on;
plot(t_n,Xe128r130_50,'r','LineWidth',1.5);
% legend('Analytical Solution with constant D&R','Mean Numerical Solution with changing D&R');
ylabel('^{128}Xe/^{130}Xe');
xlabel('Time (Gyr)');
set(gca,'FontUnits','points','FontSize',7,'FontWeight','bold',...
    'LineWidth',1.5,'LabelFontSizeMultiplier',1.0);
text(-0.126,1.01,'b','Units','normalized','FontWeight','bold','FontSize',8);
xlim([0 4.6]);

gca=subplot(2,2,2);
plot(tt,Xe132,'k','LineWidth',1.5);
hold on;
plot(t_n,Xe132_50,'r','LineWidth',1.5);
% legend('Analytical Solution with constant D&R','Mean Numerical Solution with changing D&R');
ylabel('^{132}Xe Concentration (atoms/g)');
set(gca,'FontUnits','points','FontSize',7,'FontWeight','bold',...
    'LineWidth',1.5,'LabelFontSizeMultiplier',1.0);
text(-0.126,1.01,'c','Units','normalized','FontWeight','bold','FontSize',8);
xlim([0 4.6]);

gca=subplot(2,2,4);
plot(tt,Xe132./Xe130,'k','LineWidth',1.5);
hold on;
plot(t_n,1./Xe130r132_50,'r','LineWidth',1.5);
% legend('Analytical Solution with constant D&R','Mean Numerical Solution with changing D&R','Location','southeast');
ylabel('^{132}Xe/^{130}Xe');
ylim([6 7]);
xlabel('Time (Gyr)');
set(gca,'FontUnits','points','FontSize',7,'FontWeight','bold',...
    'LineWidth',1.5,'LabelFontSizeMultiplier',1.0);
text(-0.126,1.01,'d','Units','normalized','FontWeight','bold','FontSize',8);
xlim([0 4.6]);

exportgraphics(gcf, 'test_figS1_analytical.pdf', ...
    'ContentType','image', ...   % vector/image
    'BackgroundColor','white', ... % 
    'Resolution',600); 

%%
%
figure;
set(gcf,'color','white','Units','centimeters','Position',[2 5 17.8 12]);


color_unit=[0.6, 0.8, 1.0];
color_unit = [0.3, 0.5, 0.9];
color_purple = [0.5, 0.2, 0.7];
color_unit3=[0.6 0.1 0.1];%red

color_unit1=[0.45, 0.20, 0.05];% dark red 
color_unit2=[0.8, 0.4, 0];%orange
color_unit3=[0.6 0.1 0.1];%red


% cc Formation age
dt = 0.001;% length of each timestep
t = 0:dt:tmax;
% F_5 = cc_output_margin.F_cc_5;
% F_95 = cc_output_margin.F_cc_95;
F_5 = cc_output_margin.F_Xe_5;
F_95 = cc_output_margin.F_Xe_95;
% t = t';
% t = t(:)'; 
gca=subplot(2,3,1);                       
x_fill = [t, fliplr(t)];
y_fill0 = [F_lower, fliplr(F_upper)];
y_fill1 = [F_5, fliplr(F_95)];
y_fill2 = [F_25, fliplr(F_75)];

%fill(x_fill, y_fill0, color_unit, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on;
fill(x_fill, y_fill1, color_purple, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on;
fill(x_fill, y_fill2, color_purple, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on;
p=plot(t, F_Jun_same,'-','color',[1, 0.3, 0],'LineWidth',2);
xlabel('Time (Gyr)','FontSize',6);
ylabel('CC Formation Age Distribution','FontSize',6);
% grid on;
xlim([0 4.6]);ylim([0 1.0]);
set(gca,'FontUnits','points','FontSize',8,'LabelFontSizeMultiplier',1.0,...
    'Linewidth',1,'FontWeight','bold');
text(-0.27,1.0445,'a','Units','normalized','FontWeight','bold','FontSize',10);
x=[0.2,0.22];y=[0.76,0.7];
annotation('textarrow',x,y,'String','Korenaga (2018)',...
    'HeadLength',4,'HeadWidth',4,'FontSize',6);

% h1 = patch(NaN, NaN, color_unit, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
% h2 = patch(NaN, NaN, color_unit, 'FaceAlpha', 0.75, 'EdgeColor', 'none');
% legend([p h2 h1], {'Model mean', '1σ uncertainty','2σ uncertainty'}, 'Location','northwest')


% surface age
gca=subplot(2,3,2);
x_fill = [t, fliplr(t)];
y_fill0 = [S_lower, fliplr(S_upper)];
y_fill1 = [S_5, fliplr(S_95)];
y_fill2 = [S_25, fliplr(S_75)];
%fill(x_fill, y_fill0, color_unit, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on;
fill(x_fill, y_fill1, color_purple, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on;
fill(x_fill, y_fill2, color_purple, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on;
plot(t, S_Jun_same,'-','color',[1, 0.3, 0],'LineWidth',2);
xlabel('Time (Gyr)');
ylabel('CC Surface Age Distribution');
% grid on;
xlim([0 4.6]);ylim([0 1.0]);
set(gca,'FontUnits','points','FontSize',8,'FontWeight','bold',...
    'Linewidth',1.5,'LabelFontSizeMultiplier',1.0);
text(-0.27,1.0445,'b','Units','normalized','FontWeight','bold','FontSize',10);
x=[0.49,0.51];y=[0.77,0.68];
annotation('textarrow',x,y,'String','Roberts and Spencer (2015)',...
    'HeadLength',4,'HeadWidth',4,'FontSize',6);


%
% potential temperature
default_input_for_degassing;
t=Simulation_input.time_series;
t=t/yr_s/1e9;% to Ga
t=t';
gca=subplot(2,3,3);
x_fill = [t, fliplr(t)];
y_fill0 = [Tp_lower,(fliplr(Tp_upper))];
y_fill1 = [Tp_5,(fliplr(Tp_95))];
y_fill2 = [Tp_25,(fliplr(Tp_75))];
% fill(x_fill, y_fill0,color_unit, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on;
fill(x_fill, y_fill1,color_unit, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on;
fill(x_fill, y_fill2,color_unit, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on;

% plot(t, Tp_50, 'm-', 'LineWidth', 1.5);
scatter( t_Herz,Tp_Herz, 20, ...
    'o', 'MarkerEdgeColor', [1, 0.5, 0], ...
    'MarkerFaceColor', [1, 0.3, 0], 'LineWidth', 1.0);

scatter( t_GF_fin,Tp_GF_fin, 20, ...
    'o', 'MarkerEdgeColor', [0, 0, 0], ...
    'MarkerFaceColor', [0, 0, 0], 'LineWidth', 1.0);

% % 
% stairs(t, F_50, 'Color', [1 0.4 0.2], 'LineWidth', 2); 
%
xlabel('Time (Gyr)');
ylabel('Mantle Temperature (°C)');
% grid on;
xlim([0 4.6]);
%ylim([0,1]);
% legend('accepted models','','median model','observations','location','best');
%text(2.1, 0.6, 'Korenaga (2018)', 'FontSize', 10, 'Rotation', 0);
%annotation('textarrow', [0.55 0.48], [0.55 0.65], 'String', '');
set(gca,'FontUnits','points','FontSize',8,'FontWeight','bold',...
    'Linewidth',1,'LabelFontSizeMultiplier',1.0);
text(-0.27,1.0445,'c','Units','normalized','FontWeight','bold','FontSize',10);
x=[0.81,0.814];y=[0.83,0.77];
annotation('textarrow',x,y,'String','Ganne and Feng (2017)',...
    'HeadLength',4,'HeadWidth',4,'FontSize',6);
x=[0.7663,0.7718];y=[0.6359,0.6971];
annotation('textarrow',x,y,'String','Herzberg et al. (2010)',...
    'HeadLength',4,'HeadWidth',4,'FontSize',6);


% Xe concentration
% figure;
% set(gcf,'color','white','Units','centimeters','Position',[2 5 17.8 12]);

ax_main=subplot(2,3,4);
% ax_main = subplot(2,3,4);
% ax_main.FontSize=6;
% ax_main.LabelFontSizeMultiplier = 1.0;

% set(ax_main,'FontUnits','points','FontSize',6,'LabelFontSizeMultiplier',1.0);
%set(get(gca,))

x_fill = [t, fliplr(t)];
% y_fill0 = [Xe_lower,(fliplr(Xe_upper))];
y_fill1 = [Xe_5,(fliplr(Xe_95))];
y_fill2 = [Xe_25,(fliplr(Xe_75))];

% fill(x_fill, y_fill0,color_purple, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on;
fill(x_fill, y_fill1,color_unit3, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on;
fill(x_fill, y_fill2,color_unit3, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on;

xlabel('Time (Gyr)');
ylabel('Mantle ^{130}Xe Content (atoms/g)');
% grid on;
xlim([0 4.6]);
set(ax_main,'FontUnits','points','FontSize',8,'FontWeight','bold',...
    'LineWidth',1,'LabelFontSizeMultiplier',1.0);
text(-0.27,1.0445,'d','Units','normalized','FontWeight','bold','FontSize',10);


color_unit1=[0.45, 0.20, 0.05];%dark red
tt=0:0.01:4.6;
% plot analytical solution on (need to run sup_figures.m)
hold on;
plot(tt,Xe130,'color',color_unit1,'LineWidth',1.5);
x=[0.222,0.183];y=[0.235,0.209];
annotation('textarrow',x,y,'String','Analytical Solution',...
    'HeadLength',4,'HeadWidth',4,'FontSize',6);


% zoom figure
hold(ax_main, 'on');
ax_inset = axes('Position', get(ax_main,'Position')); % 
ax_inset.Position = [ax_inset.Position(1)+0.1, ... % distance to the right
                     ax_inset.Position(2)+0.2, ... % distance to the above 
                     0.1, 0.1];                     % size-relative to the whole figure
fill(x_fill, y_fill1,color_unit3, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on;
fill(x_fill, y_fill2,color_unit3, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on;
xlim([0 4.6]);% grid on;
ylim([0 1.5e6]);
xlabel('Time (Gyr)');
% ylabel('Mantle Xe Concentration (atoms/g)');
set(ax_inset,'FontUnits','points','FontSize',6,'LabelFontSizeMultiplier',1.0);
x=[0.314,0.344];y=[0.164,0.115];
annotation('textarrow',x,y,'String','Parai and Mukhopadhyay (2018)',...
    'HeadLength',4,'HeadWidth',4,'FontSize',6);
x=[0.346,0.346];y=[0.0537,0.1121];
annotation('textarrow',x,y,'String','Present day',...
    'HeadLength',4,'HeadWidth',4,'FontSize',7);


% Xe ratio1: 128Xe/132Xe vs 130Xe/132Xe
% time=tmax-t;
time=t;



x1_atm=Xe128r132_atm(:,1);
x2_atm=Xe130r132_atm(:,1);


x1=Xe128r132_50;
x2=Xe130r132_50;

x1_5=Xe128r132_5;
x2_5=Xe130r132_5;
x1_95=Xe128r132_95;
x2_95=Xe130r132_95;
x1_10_fill=[x1_5,fliplr(x1_95)]';
x2_10_fill=[x2_5,fliplr(x2_95)]';

% alpha=1 
k_10 = boundary(x1_10_fill, x2_10_fill, 0.1);   % 0~1，


x1_25=Xe128r132_25;
x2_25=Xe130r132_25;
x1_75=Xe128r132_75;
x2_75=Xe130r132_75;
x1_50_fill=[x1_25,fliplr(x1_75)];
x2_50_fill=[x2_25,fliplr(x2_75)];



% present day mantle: a range
x1obs=Xe_obs.Xe128v132;x2obs=Xe_obs.Xe130v132;
% starting mantle: a data ponit
x1start=Xe_start_en.Xe128v132;x2start=Xe_start_en.Xe130v132;
% present day atmosphere: a data point
x1atm_obs=Xe_atm_obs.Xe128v132;x2atm_obs=Xe_atm_obs.Xe130v132;
% starting atmosphere: a data point
x1atm_start=Xe_atm_start.Xe128v132;x2atm_start=Xe_atm_start.Xe130v132;


gca=subplot(2,3,5);

% % evolutionary curve for successful params

% evolutionary curves 
surface([x1(:) x1(:)], [x2(:) x2(:)], ...
    [zeros(size(x1(:))) zeros(size(x2(:)))], ...
        [time(:) time(:)], 'EdgeColor', ...
        'interp', 'FaceColor', 'none', 'LineWidth', 2);
hold on;

hold on;
surface([x1_atm(:) x1_atm(:)], [x2_atm(:) x2_atm(:)], ...
    [zeros(size(x1_atm(:))) zeros(size(x2_atm(:)))], ...
        [time(:) time(:)], 'EdgeColor', ...
        'interp', 'FaceColor', 'none', 'LineWidth', 2);
% starting atmosphere
scatter(x1atm_start, x2atm_start, 20, 'o', 'MarkerFaceColor',[0.4 0.6 0.8], ...
        'MarkerEdgeColor','k', 'LineWidth',0.5);hold on;
% starting mantle
scatter(x1start, x2start, 20, 'o', 'MarkerFaceColor',[0.4 0.4 0.4], ...
        'MarkerEdgeColor','k', 'LineWidth',0.5);hold on;
% present day atmosphere
scatter(x1atm_obs, x2atm_obs, 20, 'o', 'MarkerFaceColor',[0.6 0.9 1], ...
        'MarkerEdgeColor','k', 'LineWidth',0.5);hold on;
% present day mantle
fill([x1obs(1) x1obs(2) x1obs(2) x1obs(1)], [x2obs(1) x2obs(1) x2obs(2) x2obs(2)], [0.5 0.5 0.5], ...
     'FaceAlpha', 0.5, 'EdgeColor', 'none');hold on;


box on;
xlim([0.060 0.085]);ylim([0.12 0.17]);


xlabel('\color[rgb]{0.6, 0.1, 0.1}^{128}Xe \color{black}/ \color[rgb]{0.9, 0.6, 0}^{132}Xe');
ylabel('\color[rgb]{0.6, 0.1, 0.1}^{130}Xe \color{black}/ \color[rgb]{0.9, 0.6, 0}^{132}Xe');
set(gca,'FontUnits','points','FontSize',8,'FontWeight','bold',...
    'Linewidth',1,'LabelFontSizeMultiplier',1.0);
text(-0.27,1.0445,'e','Units','normalized','FontWeight','bold','FontSize',10);
x=[0.5079,0.4989];y=[0.224,0.288];
annotation('textarrow',x,y,'String',['Present day mantle' newline 'Parai and Mukhopadhyay (2018)'],...
    'HeadLength',4,'HeadWidth',4,'FontSize',6);
x=[0.52,0.53];y=[0.407,0.357];
annotation('textarrow',x,y,'String','Atmosphere',...
    'HeadLength',4,'HeadWidth',4,'FontSize',6);


%
% Xe ratio2: 128Xe/132Xe vs 136Xe/132Xe

time=t;


x1_atm=Xe128r132_atm(:,1);
x2_atm=Xe136r132_atm(:,1);

x1=Xe128r132_50;
x2=Xe136r132_50;

x1_5=Xe128r132_5;
x2_5=Xe136r132_5;
x1_95=Xe128r132_95;
x2_95=Xe136r132_95;
x1_10_fill=[x1_5,fliplr(x1_95)];
x2_10_fill=[x2_5,fliplr(x2_95)];

x1_25=Xe128r132_25;
x2_25=Xe136r132_25;
x1_75=Xe128r132_75;
x2_75=Xe136r132_75;
x1_50_fill=[x1_25,fliplr(x1_75)];
x2_50_fill=[x2_25,fliplr(x2_75)];



% present day mantle: a range
x1obs=Xe_obs.Xe128v132;x2obs=Xe_obs.Xe136v132;
% starting mantle: a data ponit
x1start=Xe_start_en.Xe128v132;x2start=Xe_start_en.Xe136v132;
% present day atmosphere: a data point
x1atm_obs=Xe_atm_obs.Xe128v132;x2atm_obs=Xe_atm_obs.Xe136v132;
% starting atmosphere: a data point
x1atm_start=Xe_atm_start.Xe128v132;x2atm_start=Xe_atm_start.Xe136v132;

gca=subplot(2,3,6);
% % evolutionary curve for successful params

% evolutionary curves 
surface([x1(:) x1(:)], [x2(:) x2(:)], ...
    [zeros(size(x1(:))) zeros(size(x2(:)))], ...
        [time(:) time(:)], 'EdgeColor', ...
        'interp', 'FaceColor', 'none', 'LineWidth', 2);
hold on;


% hold on;
surface([x1_atm(:) x1_atm(:)], [x2_atm(:) x2_atm(:)], ...
    [zeros(size(x1_atm(:))) zeros(size(x2_atm(:)))], ...
        [time(:) time(:)], 'EdgeColor', ...
        'interp', 'FaceColor', 'none', 'LineWidth', 2);

% starting atmosphere
scatter(x1atm_start, x2atm_start, 20, 'o', 'MarkerFaceColor',[0.4 0.6 0.8], ...
        'MarkerEdgeColor','k', 'LineWidth',0.5);hold on;
% starting mantle
scatter(x1start, x2start, 20, 'o', 'MarkerFaceColor',[0.4 0.4 0.4], ...
        'MarkerEdgeColor','k', 'LineWidth',0.5);hold on;
% present day atmosphere
scatter(x1atm_obs, x2atm_obs, 20, 'o', 'MarkerFaceColor',[0.6 0.9 1.0], ...
        'MarkerEdgeColor','k', 'LineWidth',0.5);hold on;
% present day mantle
fill([x1obs(1) x1obs(2) x1obs(2) x1obs(1)], [x2obs(1) x2obs(1) x2obs(2) x2obs(2)], [0.5 0.5 0.5], ...
     'FaceAlpha', 0.5, 'EdgeColor', 'none');hold on;
% legend('','','','','','','modern mantle','initial mantle','modern atmosphere','initial atmosphere','Location','best');
% legend('','','','','modern mantle','initial mantle','modern atmosphere','initial atmosphere','Location','best');
h=legend('','','initial atmosphere','initial mantle','modern atmosphere','modern mantle','Location','best');
set(h,'Orientation','horizontal','Box','off',...
    'ItemTokenSize',[18 18],...
    'Position',[0.23 0.48 0.58 0.035],'FontSize',8);

box on;
xlim([0.060 0.085]);ylim([0.25 0.65]);
% colormap(flipud(jet));
colormap(glasgow);
caxis([min(time) max(time)]);
hcb=colorbar;% grid on;
title(hcb,'Gyr');
set(hcb,'Fontsize',7);
xlabel('\color[rgb]{0.6, 0.1, 0.1}^{128}Xe \color{black}/ \color[rgb]{0.9, 0.6, 0}^{132}Xe');
ylabel('\color[rgb]{0.9, 0.6, 0}^{126}Xe \color{black}/ \color[rgb]{0.9, 0.6, 0}^{132}Xe');
% ylabel('^{136}Xe/^{132}Xe');
set(gca,'FontUnits','points','FontSize',8,'FontWeight','bold',...
    'Linewidth',1,'LabelFontSizeMultiplier',1.0);
currentPosition=gca.Position;
currentPosition(3)=ax_main.Position(3);
gca.Position=currentPosition;
text(-0.27,1.0445,'f','Units','normalized','FontWeight','bold','FontSize',10);
x=[0.783,0.775];y=[0.292,0.209];
annotation('textarrow',x,y,'String',['Present day mantle' newline 'Parai and Mukhopadhyay (2018)'],...
    'HeadLength',4,'HeadWidth',4,'FontSize',6);
x=[0.796,0.815];y=[0.138,0.16];
annotation('textarrow',x,y,'String','Atmosphere',...
    'HeadLength',4,'HeadWidth',4,'FontSize',6);

% exportgraphics(gcf, 'test_fig1_analy.pdf', ...
%     'ContentType','image', ...   % vector/image
%     'BackgroundColor','white', ... % 
%     'Resolution',600); 

%% figure #3: importance & distribution of params
% clc;clear;

% % load data

%
% load data
Xe_rf_struct=load('Xe/rf_classification_results.mat');
Xe_rf_importance=Xe_rf_struct.RF_importance;
Xe_rf_target={'^{130}Xe','^{128}Xe/^{130}Xe',...
    '^{130}Xe/^{132}Xe','^{131}Xe/^{132}Xe',...
    '^{134}Xe/^{132}Xe','^{136}Xe/^{132}Xe'};
% Xe_rf_params=Xe_rf_struct.param_names;

Xe_rf_params= {'t_{s}','\kappa_g','\kappa_r','R_s','R_p',...
    'U_{BSE}','Q_c', 'T_i', '\eta_{ref}',...
    'Fr',...
    'Fd_{M}','Fd_{P}'};


% load data
success_F_struct=load('cc/out_cc_grwoth_success_lhs_MC_Krw_lowf.mat');
success_F=success_F_struct.success_F;
[size_F,nparam_cc]=size(success_F);

success_T_struct=load('thermal/out_thermal_success_lhs_test_fixK_U_MC_lowerT_Rac1100_lowf.mat');
success_T=success_T_struct.success_T;
[size_T,nparam_T]=size(success_T);

success_Xe_struct=load('Xe/out_Xe_success_lhs_new_3e7_test_fixK_U_MC_negativFD_Rac1100_new_lowf.mat');
success_Xe=success_Xe_struct.success_Xe;
[size_Xe,nparam_Xe]=size(success_Xe);




%
figure;
set(gcf,'color','white','Units','centimeters','Position',[2 5 17.8 12]);

% Feature importance of ^{130}Xe l=1
gca=subplot(3,2,1);
l=1;
b=bar(100*Xe_rf_importance(l, :)./sum(Xe_rf_importance(l, :))); 
b.FaceColor=color_unit3;
xticks(1:length(Xe_rf_params));
set(gca, 'XTickLabel', Xe_rf_params);
ylabel('Feature Importance (%)');
title(['Feature Importance for Matching ', Xe_rf_target{l}], 'FontSize', 16,'FontWeight', 'bold');
% grid on;
set(gca,'FontUnits','points','FontSize',7,'FontWeight','bold',...
    'Linewidth',1,'LabelFontSizeMultiplier',1.0);
set(gca,'Position',[0.06 0.7177 0.505 0.2]);
text(-0.1,1.01,'a','Units','normalized','FontWeight','bold','FontSize',8);

% Feature importance of ^{128}Xe/^{130}Xe l=2
gca=subplot(3,2,3);
l=2;
b=bar(100*Xe_rf_importance(l, :)./sum(Xe_rf_importance(l, :))); 
b.FaceColor=color_unit3;
xticks(1:length(Xe_rf_params));
set(gca, 'XTickLabel', Xe_rf_params);
ylabel('Feature Importance (%)');
title(['Feature Importance for Matching ', Xe_rf_target{l}], 'FontSize', 16,'FontWeight', 'bold');
% grid on;
set(gca,'FontUnits','points','FontSize',7,'FontWeight','bold',...
    'Linewidth',1,'LabelFontSizeMultiplier',1.0);
set(gca,'Position',[0.06 0.4181 0.505 0.2]);
text(-0.1,1.01,'b','Units','normalized','FontWeight','bold','FontSize',8);

% Feature importance of '^{130}Xe/^{132}Xe' l=3
gca=subplot(3,2,5);
l=6;
b=bar(100*Xe_rf_importance(l, :)./sum(Xe_rf_importance(l, :))); 
b.FaceColor=color_unit3;
xticks(1:length(Xe_rf_params));
set(gca, 'XTickLabel', Xe_rf_params);
ylabel('Feature Importance (%)');
title(['Feature Importance for Matching ', Xe_rf_target{l}], 'FontSize', 16,'FontWeight', 'bold');
% grid on;
set(gca,'FontUnits','points','FontSize',7,'FontWeight','bold',...
    'Linewidth',1,'LabelFontSizeMultiplier',1.0);
set(gca,'Position',[0.06 0.1184 0.505 0.2]);
text(-0.1,1.01,'c','Units','normalized','FontWeight','bold','FontSize',8);


% histogram of Ubse
k=6;
y2=success_T(:,k)*1e9;
y3=success_Xe(:,k)*1e9;

gca=subplot(3,2,2);
histogram(y2, 15, 'Normalization','probability',...            
              'EdgeColor',color_unit, ...      % boundary
              'FaceColor','none', ...   % fill
              'LineWidth',1.2);         % boundary
hold on;% 'BinWidth',0.1, ...
histogram(y3, 15, 'Normalization','probability',...
              'EdgeColor',color_unit3, ...      % boundary
              'FaceColor','none', ...   % fill
              'LineWidth',1.2);         % boundary
ylim([0 0.6]);
yt = get(gca,'YTick');            % 
set(gca, 'YTick', yt, 'YTickLabel', yt*100);
set(gca,'Position',[0.613 0.717 0.351 0.198]);

xlabel('U_{BSE} (ppb)');ylabel('Frequency (%)');
title('Histogram of Solutions Satisfying All Constraints', 'FontSize', 12,'FontWeight', 'bold');
set(gca,'FontUnits','points','FontSize',7,...
    'Linewidth',1,'FontWeight','bold','LabelFontSizeMultiplier',1.0);
text(-0.1,1.01,'d','Units','normalized','FontWeight','bold','FontSize',8);



% histogram of kappa_gcc or kappa_rcc
k=2;
y1=success_F(:,k);
y2=success_T(:,k);
y3=success_Xe(:,k);

gca=subplot(3,2,4);
histogram(y1, 15, 'Normalization','probability',...            
              'EdgeColor',color_purple, ...      % boundary
              'FaceColor','none', ...   % fill
              'LineWidth',1.2);         % boundary
hold on;% 'BinWidth',0.1, ...
histogram(y2, 15, 'Normalization','probability',...
              'EdgeColor',color_unit, ...      % boundary
              'FaceColor','none', ...   % fill
              'LineWidth',1.2);         % boundary
histogram(y3, 15, 'Normalization','probability',...
              'EdgeColor',color_unit3, ...      % boundary
              'FaceColor','none', ...   % fill
              'LineWidth',1.2);         % boundary
ylim([0 0.6]);
yt = get(gca,'YTick');            % 
set(gca, 'YTick', yt, 'YTickLabel', yt*100);
set(gca,'Position',[0.613 0.4181 0.351 0.198]);


% legend('Stage 1','Stage 2','Stage 3','FontSize',6);
xlabel('\kappa_g (Gyr^{-1})');ylabel('Frequency (%)');
set(gca,'FontUnits','points','FontSize',7,...
    'Linewidth',1,'FontWeight','bold','LabelFontSizeMultiplier',1.0);
text(-0.1,1.01,'e','Units','normalized','FontWeight','bold','FontSize',8);

% % for kappa_gcc: zoom figure



% histogram of Rs
k=4;
y1=success_F(:,k)/1e22;
y2=success_T(:,k)/1e22;
y3=success_Xe(:,k)/1e22;

gca=subplot(3,2,6);
histogram(y1, 15, 'Normalization','probability',...            
              'EdgeColor',color_purple, ...      % boundary
              'FaceColor','none', ...   % fill
              'LineWidth',1.2);         % boundary
hold on;% 'BinWidth',0.1, ...
histogram(y2, 15, 'Normalization','probability',...
              'EdgeColor',color_unit, ...      % boundary
              'FaceColor','none', ...   % fill
              'LineWidth',1.2);         % boundary
histogram(y3, 15, 'Normalization','probability',...
              'EdgeColor',color_unit3, ...      % boundary
              'FaceColor','none', ...   % fill
              'LineWidth',1.2);         % boundary

ylim([0 0.6]);
yt = get(gca,'YTick');            % 
set(gca, 'YTick', yt, 'YTickLabel', yt*100);
set(gca,'Position',[0.613 0.1184 0.351 0.198]);


           % 
xlabel('R_{s} (e22 kg/Gyr)');ylabel('Frequency (%)');
set(gca,'FontUnits','points','FontSize',7,'FontWeight','bold',...
    'Linewidth',1,'LabelFontSizeMultiplier',1.0);
text(-0.1,1.01,'f','Units','normalized','FontWeight','bold','FontSize',8);


h=legend('CC','Thermal','Xe');
set(h,'Orientation','horizontal','Box','off',...
    'ItemTokenSize',[18 18],...
    'Position',[0.645 0.87 0.58 0.035],'FontSize',6);

% exportgraphics(gcf, 'test_fig2.pdf', ...
%     'ContentType','image', ...   % vector/image
%     'BackgroundColor','white', ... % 
%     'Resolution',600); 

%% figure #2: cc curves after vs before; degassing analysis
% clc;clear;
par_solver_new;
tmax = t_pd/yr_s/1e9;
dt = 0.001;% length of each timestep
t = 0:dt:tmax;
nt=length(t);

% load data

cc_output_margin_struct=load('Xe/cc_output_for_figure_margin_fin0_lowf_all.mat');
cc_output_margin=cc_output_margin_struct.cc_output_margin;


Xe_output_margin_struct=load('Xe/Xe_output_for_figure_margin_new_lowf.mat');
Xe_output_margin=Xe_output_margin_struct.Xe_output_margin;


Xe_output_ex_struct=load('Xe/Xe_output_for_figure_ex_new_lowf.mat');
Xe_output_ex=Xe_output_ex_struct.Xe_output;

% 90% successful

NGcc_cc_5 = cc_output_margin.NGcc_cc_5;
NGcc_cc_95 = cc_output_margin.NGcc_cc_95;


Gcc_cc_5 = cc_output_margin.Gcc_cc_5;
Gcc_cc_95 = cc_output_margin.Gcc_cc_95;

Rcc_cc_5 = cc_output_margin.Rcc_cc_5;
Rcc_cc_95 = cc_output_margin.Rcc_cc_95;

cc_output_margin_struct=load('Xe/cc_output_for_figure_margin_fin0_lowf_ccRcc.mat');
cc_output_margin=cc_output_margin_struct.cc_output_margin;
Rcc_cc_lower = cc_output_margin.Rcc_lower;
Rcc_cc_upper = cc_output_margin.Rcc_upper;

cc_output_margin_struct=load('Xe/cc_output_for_figure_margin_fin0_lowf.mat');
cc_output_margin=cc_output_margin_struct.cc_output_margin;


NGcc_Xe_5 = cc_output_margin.NGcc0_Xe_5;
NGcc_Xe_95 = cc_output_margin.NGcc0_Xe_95;


Gcc_Xe_5 = cc_output_margin.Gcc0_Xe_5;
Gcc_Xe_95 = cc_output_margin.Gcc0_Xe_95;


Rcc_Xe_5 = cc_output_margin.Rcc0_Xe_5;
Rcc_Xe_95 = cc_output_margin.Rcc0_Xe_95;

%

Dcc_5 = Xe_output_margin.D0_cc_5;
Dcc_50 = Xe_output_margin.D0_cc_50;
Dcc_95 = Xe_output_margin.D0_cc_95;


Dmor_5 = Xe_output_margin.D0_mor_5;
Dmor_50 = Xe_output_margin.D0_mor_50;
Dmor_95 = Xe_output_margin.D0_mor_95;


Dplume_5 = Xe_output_margin.D0_p_5;
Dplume_50 = Xe_output_margin.D0_p_50;
Dplume_95 = Xe_output_margin.D0_p_95;


%
figure;
set(gcf,'color','white','Units','centimeters','Position',[2 5 11.4 12]);
color_unit=[0.6, 0.8, 1.0];
color_purple = [0.5, 0.2, 0.7];

% % NGcc
gca=subplot(2,2,1);
% text(0.02,0.98,'(a)','Units','normalized','FontWeight','bold','FontSize',8);
x_fill = [t, fliplr(t)];


y_fill1 = [NGcc_Xe_5,(fliplr(NGcc_Xe_95))];


z_fill1 = [NGcc_cc_5,(fliplr(NGcc_cc_95))];

 fill(x_fill, z_fill1,color_purple, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on;

 fill(x_fill, y_fill1,color_unit3, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on;

xlabel('Time (Gyr)'); ylabel('Net CC Growth (kg)'); % grid on;
xlim([0 4.6]);

h1 = patch(NaN, NaN, color_purple, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
h2 = patch(NaN, NaN, [0.567, 0.133, 0.300], 'FaceAlpha', 0.75, 'EdgeColor', 'none');
h=legend([h1 h2], {'CC','Xe'}, 'Location','southeast');
set(gca,'FontUnits','points','FontSize',7,'FontWeight','bold',...
    'Linewidth',1,'LabelFontSizeMultiplier',1.0);

text(-0.18,1.01,'a','Units','normalized','FontWeight','bold','FontSize',8);
x=[0.233,0.168];y=[0.715,0.812];

% % Rcc
gca=subplot(2,2,2);
x_fill = [t, fliplr(t)];

y_fill1 = [Rcc_Xe_5,(fliplr(Rcc_Xe_95))];

z_fill1 = [Rcc_cc_5,(fliplr(Rcc_cc_95))];

 fill(x_fill, z_fill1,color_purple, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on;

 fill(x_fill, y_fill1,color_unit3, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on;
xlabel('Time (Gyr)'); ylabel('CC Recycling Rate (kg/Gyr)'); % grid on;
xlim([0 4.6]);
set(gca,'FontUnits','points','FontSize',7,'FontWeight','bold',...
    'Linewidth',1,'LabelFontSizeMultiplier',1.0);
text(-0.18,1.01,'b','Units','normalized','FontWeight','bold','FontSize',8);

% % Gcc
gca=subplot(2,2,3);
x_fill = [t, fliplr(t)];

y_fill1 = [Gcc_Xe_5,(fliplr(Gcc_Xe_95))];

z_fill1 = [Gcc_cc_5,(fliplr(Gcc_cc_95))];

 fill(x_fill, z_fill1,color_purple, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on;

 fill(x_fill, y_fill1,color_unit3, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on;
xlabel('Time (Gyr)'); ylabel('CC Generation Rate (kg/Gyr)'); % grid on;
xlim([0 4.6]);
set(gca,'FontUnits','points','FontSize',7,'FontWeight','bold',...
    'Linewidth',1,'LabelFontSizeMultiplier',1.0);
text(-0.18,1.01,'c','Units','normalized','FontWeight','bold','FontSize',8);



% % Degassing analysis ()
default_input_for_degassing;
t=Simulation_input.time_series;
t=t/yr_s/1e9;% to Ga
t=t';

gca=subplot(2,2,4);

[maxc,maxc_loc]=max(Dcc_50);

semilogy(t(maxc_loc:end),Dcc_50(maxc_loc:end)*yr_s,'color',color_purple,'Linewidth',1);hold on;
semilogy(t,Dmor_50*yr_s,'b','Linewidth',1);hold on;
semilogy(t,Dplume_50*yr_s,'color',color_unit,'Linewidth',1);hold on;

xlabel('Time (Gyr)'); ylabel('Degassing flux (kg/yr)'); % grid on;
xlim([0 4.6]);
ylim([1e14 3e16]);
h=legend('D_{CC}','D_{MOR}','D_{Plume}','Location','northeast');
set(h,...
    'ItemTokenSize',[12 12],...
    'FontSize',6);
set(gca,'FontUnits','points','FontSize',7,'FontWeight','bold',...
    'Linewidth',1.5,'LabelFontSizeMultiplier',1.0);
text(-0.18,1.01,'d','Units','normalized','FontWeight','bold','FontSize',8);

% exportgraphics(gcf, 'test_fig3.pdf', ...
%     'ContentType','image', ...   % vector/image
%     'BackgroundColor','white', ... % 
%     'Resolution',600); 


%% % figure #4: why highcc and lowU
%% highcc 
% clc;clear all;
% contour_struct=load('Xe/fin_new_2000_RpFr_Xe_3e7.mat');
contour_struct=load('Xe/fin_new_2007_RpFr_Xe_3e7.mat');

Ubse_Fr_Xe_contours=contour_struct.Ubse_Fr_Xe_contours;
Rs_Fr_Xe_contours=contour_struct.Rs_Fr_Xe_contours;
kgcc_Fr_Xe_contours=contour_struct.kgcc_Fr_Xe_contours;
col_names=contour_struct.col_names;
Xq=contour_struct.FrXe;

figure;
set(gcf,'color','white','Units','centimeters','Position',[2 5 11.4 12]);

gca=subplot(3,2,1);
% histogram of kappa_gcc or kappa_rcc
k=2;
y1=success_F(:,k);
y2=success_T(:,k);
y3=success_Xe(:,k);

histogram(y1, 15, 'Normalization','probability',...            
              'EdgeColor',color_purple, ...      % boundary
              'FaceColor','none', ...   % fill
              'LineWidth',1.2);         % boundary
hold on;% 'BinWidth',0.1, ...
histogram(y2, 15, 'Normalization','probability',...
              'EdgeColor',color_unit, ...      % boundary
              'FaceColor','none', ...   % fill
              'LineWidth',1.2);         % boundary
histogram(y3, 15, 'Normalization','probability',...
              'EdgeColor',color_unit3, ...      % boundary
              'FaceColor','none', ...   % fill
              'LineWidth',1.2);         % boundary
ylim([0 0.6]);
yt = get(gca,'YTick');            % 
set(gca, 'YTick', yt, 'YTickLabel', yt*100);
% set(gca,'Position',[0.613 0.4181 0.351 0.198]);


legend('CC','Thermal','Xe','FontSize',6);
xlabel('\kappa_g (Gyr^{-1})');ylabel('Frequency (%)');
set(gca,'FontUnits','points','FontSize',7,...
    'Linewidth',1,'FontWeight','bold','LabelFontSizeMultiplier',1.0);
text(-0.2,1.01,'a','Units','normalized','FontWeight','bold','FontSize',8);


gca=subplot(3,2,3);
% histogram of Ubse
k=6;
y2=success_T(:,k)*1e9;
y3=success_Xe(:,k)*1e9;

histogram(y2, 15, 'Normalization','probability',...            
              'EdgeColor',color_unit, ...      % boundary
              'FaceColor','none', ...   % fill
              'LineWidth',1.2);         % boundary
hold on;% 'BinWidth',0.1, ...
histogram(y3, 15, 'Normalization','probability',...
              'EdgeColor',color_unit3, ...      % boundary
              'FaceColor','none', ...   % fill
              'LineWidth',1.2);         % boundary
ylim([0 0.6]);
yt = get(gca,'YTick');            % 
set(gca, 'YTick', yt, 'YTickLabel', yt*100);
% set(gca,'Position',[0.613 0.717 0.351 0.198]);

xlabel('U_{BSE} (ppb)');ylabel('Frequency (%)');
% title('Histogram of Solutions Satisfying All Constraints', 'FontSize', 12,'FontWeight', 'bold');
set(gca,'FontUnits','points','FontSize',7,...
    'Linewidth',1,'FontWeight','bold','LabelFontSizeMultiplier',1.0);
text(-0.2,1.01,'c','Units','normalized','FontWeight','bold','FontSize',8);




gca=subplot(3,2,2);
ax = gca;
%
results=kgcc_Fr_Xe_contours;
% results=Rs_Fr_Xe_contours;
Yq=contour_struct.kappa_gcc;
Yq= -1:3.2:31;
% Yq= contour_struct.Rs;
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

color_unit1=[0.45, 0.20, 0.05];%dark red
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
% subplot(2,2,1);

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
xlabel('Fr', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('\kappa_g', 'FontSize', 12, 'FontWeight', 'bold')
% title('test',...
%     'FontSize', 14, 'FontWeight', 'bold')%Xeini=3.2e7atoms/gram, 
set(gca,'FontUnits','points','FontSize',7,'FontWeight','bold',...
    'Linewidth',1.5,'LabelFontSizeMultiplier',1.0);
box on;
text(-0.2,1.01,'b','Units','normalized','FontWeight','bold','FontSize',8);



% low U 
% 

col_names=contour_struct.col_names;
Xq=contour_struct.FrXe;


%
gca=subplot(3,2,4);
ax = gca;
results=Ubse_Fr_Xe_contours;
Yq=contour_struct.Ubse;
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

color_unit1=[0.45, 0.20, 0.05];% dark red 
color_unit2=[0.8, 0.4, 0];%orange
color_unit3=[0.6 0.1 0.1];%red

hPatch1 = patch(NaN,NaN,color_unit1,'FaceAlpha',alpha*2,'EdgeColor', color_unit1, 'LineWidth', 1.5, 'LineStyle', '--');
hPatch2 = patch(NaN,NaN,color_unit3,'FaceAlpha',alpha*2,'EdgeColor', color_unit3, 'LineWidth', 1.5, 'LineStyle', '--');
hPatch3 = patch(NaN,NaN,color_unit2,'FaceAlpha',alpha*2,'EdgeColor','none');


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
% ax = axes;
alpha=0.2;
minX=min(min(Xq));maxX=max(max(Xq));
minY=min(min(Yq));maxY=max(max(Yq));


himg1=image(ax, img1, 'AlphaData', alpha*2,'XData', [minX,maxX], 'YData', [maxY,minY]); % 设置透明度
hold on;
himg2=image(ax, img2, 'AlphaData', alpha,'XData', [minX,maxX], 'YData', [maxY,minY]);
axis(ax, 'tight');
hold on;
himg3=image(ax, img3, 'AlphaData', alpha,'XData', [minX,maxX], 'YData', [maxY,minY]);
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
xlabel('Fr', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('U_{BSE}', 'FontSize', 12, 'FontWeight', 'bold')
% title('test',...
%     'FontSize', 14, 'FontWeight', 'bold') %Xeini=3.2e7atoms/gram, 
set(gca,'FontUnits','points','FontSize',7,'FontWeight','bold',...
    'Linewidth',1,'LabelFontSizeMultiplier',1.0);
h=legend('^{130}Xe concentration','^{128}Xe/^{130}Xe','^{131,132,134,136}Xe/^{130}Xe');
% set(h,'Orientation','horizontal','Box','off',...
%     'ItemTokenSize',[18 18],...
%     'Position',[0.2 0.48 0.58 0.035],'FontSize',6);
set(h,'Orientation','horizontal','Box','off',...
    'ItemTokenSize',[18 18],...
    'Position',[0.1309 0.9541 0.75 0.04],'FontSize',6);
box on;
text(-0.2,1.01,'d','Units','normalized','FontWeight','bold','FontSize',8);



% ex
Xe_output_ex_struct=load('Xe/Xe_output_for_figure_ex_new_lowf.mat');
Xe_output_ex=Xe_output_ex_struct.Xe_output;

Xe128r132_atm=Xe_output_ex.Xe128r132_atm(:,1);
Xe130r132_atm=Xe_output_ex.Xe130r132_atm(:,1);
Xe136r132_atm=Xe_output_ex.Xe136r132_atm(:,1);

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

default_input_for_degassing;
t=Simulation_input.time_series;
t=t/yr_s/1e9;% to Ga
t=t';
% time=tmax-t;
time=t;

gca=subplot(3,2,5);

x1_atm=Xe128r132_atm;
x2_atm=Xe130r132_atm;

% present day mantle: a range
x1obs=Xe_obs.Xe128v132;x2obs=Xe_obs.Xe130v132;
% starting mantle: a data ponit
x1start=Xe_start_en.Xe128v132;x2start=Xe_start_en.Xe130v132;
% present day atmosphere: a data point
x1atm_obs=Xe_atm_obs.Xe128v132;x2atm_obs=Xe_atm_obs.Xe130v132;
% starting atmosphere: a data point
x1atm_start=Xe_atm_start.Xe128v132;x2atm_start=Xe_atm_start.Xe130v132;


x1=Xe128t./Xe132t;
x2=Xe130t./Xe132t;


x1_lowcc=Xe128t_lowcc_Rs./Xe132t_lowcc_Rs;
x2_lowcc=Xe130t_lowcc_Rs./Xe132t_lowcc_Rs;

x1_highU=Xe128t_highU./Xe132t_highU;
x2_highU=Xe130t_highU./Xe132t_highU;


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
surface([x1_atm(:) x1_atm(:)], [x2_atm(:) x2_atm(:)], ...
    [zeros(size(x1_atm(:))) zeros(size(x2_atm(:)))], ...
        [time(:) time(:)], 'EdgeColor', ...
        'interp', 'FaceColor', 'none', 'LineWidth', 2);

% starting atmosphere
scatter(x1atm_start, x2atm_start, 20, 'o', 'MarkerFaceColor',[0.4 0.6 0.8], ...
        'MarkerEdgeColor','k', 'LineWidth',0.5);hold on;
% starting mantle
scatter(x1start, x2start, 20, 'o', 'MarkerFaceColor',[0.4 0.4 0.4], ...
        'MarkerEdgeColor','k', 'LineWidth',0.5);hold on;
% present day atmosphere
scatter(x1atm_obs, x2atm_obs, 20, 'o', 'MarkerFaceColor',[0.6 0.9 1.0], ...
        'MarkerEdgeColor','k', 'LineWidth',0.5);hold on;
% present day mantle
fill([x1obs(1) x1obs(2) x1obs(2) x1obs(1)], [x2obs(1) x2obs(1) x2obs(2) x2obs(2)], [0.5 0.5 0.5], ...
     'FaceAlpha', 0.5, 'EdgeColor', 'none');hold on;
% legend('','','','','initial atmosphere','initial mantle','modern atmosphere','modern mantle','Location','best');
box on;

xlim([0.060 0.085]);ylim([0.12 0.17]);
colormap(glasgow);
% hcb=colorbar;% grid on;
% title(hcb,'Gyr');
xlabel('^{128}Xe/^{132}Xe');
ylabel('^{130}Xe/^{132}Xe');
set(gca,'FontUnits','points','FontSize',7,'FontWeight','bold',...
    'Linewidth',1,'LabelFontSizeMultiplier',1.0);
gca_n=gca;
text(-0.2,1.01,'e','Units','normalized','FontWeight','bold','FontSize',8);
% x=[0.344,0.307];y=[0.218,0.294];
x=[0.3313,0.3063];y=[0.2029,0.2353];
annotation('textarrow',x,y,'String','Slow CC growth',...
    'HeadLength',4,'HeadWidth',4,'FontSize',6);
% x=[0.229,0.272];y=[0.209,0.253];
x=[0.229,0.272];y=[0.1737,0.2177];
annotation('textarrow',x,y,'String','High U_{BSE}',...
    'HeadLength',4,'HeadWidth',4,'FontSize',6);
% hcb5 = colorbar; 
% hcb5.Visible = 'off';

%
gca=subplot(3,2,6);

x1_atm=Xe128r132_atm;
x2_atm=Xe136r132_atm;

% present day mantle: a range
x1obs=Xe_obs.Xe128v132;x2obs=Xe_obs.Xe136v132;
% starting mantle: a data ponit
x1start=Xe_start_en.Xe128v132;x2start=Xe_start_en.Xe136v132;
% present day atmosphere: a data point
x1atm_obs=Xe_atm_obs.Xe128v132;x2atm_obs=Xe_atm_obs.Xe136v132;
% starting atmosphere: a data point
x1atm_start=Xe_atm_start.Xe128v132;x2atm_start=Xe_atm_start.Xe136v132;


x1=Xe128t./Xe132t;
x2=Xe136t./Xe132t;


x1_lowcc=Xe128t_lowcc_Rs./Xe132t_lowcc_Rs;
x2_lowcc=Xe136t_lowcc_Rs./Xe132t_lowcc_Rs;

x1_highU=Xe128t_highU./Xe132t_highU;
x2_highU=Xe136t_highU./Xe132t_highU;


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

hold on;
surface([x1_atm(:) x1_atm(:)], [x2_atm(:) x2_atm(:)], ...
    [zeros(size(x1_atm(:))) zeros(size(x2_atm(:)))], ...
        [time(:) time(:)], 'EdgeColor', ...
        'interp', 'FaceColor', 'none', 'LineWidth', 2);

% starting atmosphere
scatter(x1atm_start, x2atm_start, 20, 'o', 'MarkerFaceColor',[0.4 0.6 0.8], ...
        'MarkerEdgeColor','k', 'LineWidth',0.5);hold on;
% starting mantle
scatter(x1start, x2start, 20, 'o', 'MarkerFaceColor',[0.4 0.4 0.4], ...
        'MarkerEdgeColor','k', 'LineWidth',0.5);hold on;
% present day atmosphere
scatter(x1atm_obs, x2atm_obs, 20, 'o', 'MarkerFaceColor',[0.6 0.9 1.0], ...
        'MarkerEdgeColor','k', 'LineWidth',0.5);hold on;
% present day mantle
fill([x1obs(1) x1obs(2) x1obs(2) x1obs(1)], [x2obs(1) x2obs(1) x2obs(2) x2obs(2)], [0.5 0.5 0.5], ...
     'FaceAlpha', 0.5, 'EdgeColor', 'none');hold on;
h=legend('','','','','initial atmosphere','initial mantle','modern atmosphere','modern mantle','Location','best');
set(h,'Orientation','horizontal','Box','off',...
    'ItemTokenSize',[18 18],...
    'Position',[0.03 0.015 0.9 0.035],'FontSize',6);
box on;

xlim([0.060 0.085]);ylim([0.25 0.65]);
colormap(glasgow);
hcb=colorbar;% grid on;
set(hcb, 'Position', [0.93 0.1309 0.0226 0.1956]);
title(hcb,'Gyr');
% set(hcb,'Position',[1 0.115 0.027 0.33],'Fontsize',6);
xlabel('^{128}Xe/^{132}Xe');
ylabel('^{136}Xe/^{132}Xe');
set(gca,'FontUnits','points','FontSize',7,'FontWeight','bold',...
    'Linewidth',1,'LabelFontSizeMultiplier',1.0);

text(-0.2,1.01,'f','Units','normalized','FontWeight','bold','FontSize',8);
% x=[0.799,0.746];y=[0.356,0.265];
x=[0.774,0.7522];y=[0.2529,0.2091];
annotation('textarrow',x,y,'String','Slow CC growth',...
    'HeadLength',4,'HeadWidth',4,'FontSize',6);
% x=[0.665,0.703];y=[0.377,0.309];
x=[0.6756,0.7034];y=[0.2565, 0.2236];
annotation('textarrow',x,y,'String','High U_{BSE}',...
    'HeadLength',4,'HeadWidth',4,'FontSize',6);


% 
% 
% exportgraphics(gcf, 'test_fig4_new.pdf', ...
%     'ContentType','image', ...   % vector/image
%     'BackgroundColor','white', ... % 
%     'Resolution',600); 





