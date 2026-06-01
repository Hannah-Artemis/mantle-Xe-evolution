%% % sup figure2-3
% rf_cc_T; histogram_Xe; SHAP

%% S2
% load data
cc_rf_struct=load('cc/rf_classification_results_cc.mat');
cc_rf_importance_f=cc_rf_struct.RF_importance(1,:);
cc_rf_importance_s=cc_rf_struct.RF_importance(2,:);
cc_rf_params= cc_rf_struct.param_names;

thermal_rf_struct=load('thermal/rf_classification_results_T.mat');
thermal_rf_importance=thermal_rf_struct.RF_importance;
thermal_rf_target= {'Present day Tp','Historical Tp','Present day Qs'};
thermal_rf_params= {'ts','kappa_g','kappa_r','Rs','Rp','U_{BSE}','Qcpd','Ti','\eta_{ref}'};


%
figure;
set(gcf,'color','white','Units','centimeters','Position',[2 5 17.8 12]);

% Feature importance of ^{130}Xe l=1
gca=subplot(3,2,3);
b=bar(100*abs(cc_rf_importance_f)./sum(abs(cc_rf_importance_f))); 
b.FaceColor='k';
xticks(1:length(cc_rf_params));
set(gca, 'XTickLabel', cc_rf_params);
ylabel('Feature Importance(%)');
title('Influence on CC formation age', 'FontSize', 16,'FontWeight', 'bold');
% grid on;
set(gca,'FontUnits','points','FontSize',7,'FontWeight','bold',...
    'Linewidth',1.5,'LabelFontSizeMultiplier',1.0);
% set(gca,'Position',[0.07 0.7177 0.4 0.2]);
text(-0.126,1.01,'a','Units','normalized','FontWeight','bold','FontSize',8);


gca=subplot(3,2,5);
b=bar(100*cc_rf_importance_s./sum(cc_rf_importance_s)); 
b.FaceColor='k';
xticks(1:length(cc_rf_params));
set(gca, 'XTickLabel', cc_rf_params);
ylabel('Feature Importance(%)');
title('Influence on CC surface age', 'FontSize', 16,'FontWeight', 'bold');
% grid on;
set(gca,'FontUnits','points','FontSize',7,'FontWeight','bold',...
    'Linewidth',1.5,'LabelFontSizeMultiplier',1.0);
% set(gca,'Position',[0.07 0.7177 0.4 0.2]);
text(-0.126,1.01,'b','Units','normalized','FontWeight','bold','FontSize',8);


gca=subplot(3,2,2);
l=1;
b=bar(100*thermal_rf_importance(l, :)./sum(thermal_rf_importance(l, :))); 
b.FaceColor='k';
xticks(1:length(thermal_rf_params));
set(gca, 'XTickLabel', thermal_rf_params);
ylabel('Feature Importance(%)');
title(['Influence on ', thermal_rf_target{l}], 'FontSize', 16,'FontWeight', 'bold');
% grid on;
set(gca,'FontUnits','points','FontSize',7,'FontWeight','bold',...
    'Linewidth',1.5,'LabelFontSizeMultiplier',1.0);
set(gca,'Position',[0.5456 0.7177 0.4 0.2]);
text(-0.126,1.01,'c','Units','normalized','FontWeight','bold','FontSize',8);

gca=subplot(3,2,4);
l=2;
b=bar(100*thermal_rf_importance(l, :)./sum(thermal_rf_importance(l, :))); 
b.FaceColor='k';
xticks(1:length(thermal_rf_params));
set(gca, 'XTickLabel', thermal_rf_params);
ylabel('Feature Importance(%)');
title(['Influence on ', thermal_rf_target{l}], 'FontSize', 16,'FontWeight', 'bold');
% grid on;
set(gca,'FontUnits','points','FontSize',7,'FontWeight','bold',...
    'Linewidth',1.5,'LabelFontSizeMultiplier',1.0);
set(gca,'Position',[0.5456 0.4154 0.4 0.2]);
text(-0.126,1.01,'d','Units','normalized','FontWeight','bold','FontSize',8);

gca=subplot(3,2,6);
l=3;
b=bar(100*thermal_rf_importance(l, :)./sum(thermal_rf_importance(l, :))); 
b.FaceColor='k';
xticks(1:length(thermal_rf_params));
set(gca, 'XTickLabel', thermal_rf_params);
ylabel('Feature Importance(%)');
title(['Influence on ', thermal_rf_target{l}], 'FontSize', 16,'FontWeight', 'bold');
% grid on;
set(gca,'FontUnits','points','FontSize',7,'FontWeight','bold',...
    'Linewidth',1.5,'LabelFontSizeMultiplier',1.0);
set(gca,'Position',[0.5456 0.11 0.4 0.2]);
text(-0.126,1.01,'e','Units','normalized','FontWeight','bold','FontSize',8);

exportgraphics(gcf, 'test_figS2_rf.pdf', ...
    'ContentType','image', ...   % vector/image
    'BackgroundColor','white', ... % 
    'Resolution',600); 


%% % S3 - histogram
% load data
success_F_struct=load('cc/out_cc_grwoth_success_lhs_MC_Krw_lowf.mat');
success_F=success_F_struct.success_F;
[size_F,nparam_cc]=size(success_F);

success_T_struct=load('thermal/out_thermal_success_lhs_test_fixK_U_MC_lowerT_Rac1100_lowf.mat');
success_T=success_T_struct.success_T;
[size_T,nparam_T]=size(success_T);

results_Xe_struct=load('Xe/out_Xe_lhs_new_3e7_test_fixK_U_MC_lowerT_alot_Rac1100_new_lowf.mat');
results=results_Xe_struct.results;

col_names = {'ts_cc','kappa_gcc','kappa_rcc','Rs','Rp',...
    'Ubse','Qc_pd', 'dQc', 'Ti', 'eta_ref',...
    'Xm_init','Ms_init','Fr_w',...
    'Xe_init','Fr_Xe',...
    'Fd_mor','Fd_p',...
    'Er_Pmor_pd',... % 'RMSE_mor_processing rate',...
    'Er_xec','Er_xe128v130','Er_xe128v132','Er_xe130v132','Er_xe131v132','Er_xe134v132','Er_xe136v132'};

param_names = {'ts','\kappa_g','\kappa_r','Rs','Rp',...
    'Ubse','Qc_pd', 'Ti', 'eta_ref',...
    'Fr_Xe',...
    'Fd_mor','Fd_p'};

Erxec_idx = find(strcmp(col_names, 'Er_xec'));
Erxe128v130_idx = find(strcmp(col_names, 'Er_xe128v130'));
Erxe128v132_idx = find(strcmp(col_names, 'Er_xe128v132'));
Erxe130v132_idx = find(strcmp(col_names, 'Er_xe130v132'));
Erxe131v132_idx = find(strcmp(col_names, 'Er_xe131v132'));
Erxe134v132_idx = find(strcmp(col_names, 'Er_xe134v132'));
Erxe136v132_idx = find(strcmp(col_names, 'Er_xe136v132'));

resid_col_names = {'^{130}Xe','^{128}Xe/^{130}Xe','^{130}Xe/^{132}Xe','^{131}Xe/^{132}Xe','^{134}Xe/^{132}Xe','^{136}Xe/^{132}Xe'};
resid_col_ids = [Erxec_idx,Erxe128v130_idx, Erxe130v132_idx, Erxe131v132_idx, Erxe134v132_idx, Erxe136v132_idx];

% success_mask = (    abs( results(:, Erxec_idx) ) <= crit_misfit_xec & ...
%                     abs( results(:, Erxe128v130_idx) ) <= crit_misfit_xe128v130 & ... % results(:, Erxe128v132_idx) <= crit_misfit_xe128v132 & ...                    
%                     abs( results(:, Erxe130v132_idx) ) <= crit_misfit_xe130v132 & ...
%                     abs( results(:, Erxe131v132_idx) ) <= crit_misfit_xe131v132 & ...
%                     abs( results(:, Erxe134v132_idx) ) <= crit_misfit_xe134v132 & ...
%                     abs( results(:, Erxe136v132_idx) ) <= crit_misfit_xe136v132 & ...
%                     abs( results(:, ErPmor_idx) ) <= crit_misfit_Pmor );


figure;
set(gcf,'color','white','Units','centimeters','Position',[2 5 17.8 12]);

% for Xe concentration
success_mask = abs( results(:, Erxec_idx) ) <= 1;
success_Xe=results(success_mask,:);
gca=subplot(3,2,1);
% histogram of Ubse
k=6;
y2=success_T(:,k)*1e9;
y3=success_Xe(:,k)*1e9;


histogram(y2, 15, 'Normalization','probability',...            
              'EdgeColor','b', ...      % boundary
              'FaceColor','none', ...   % fill
              'LineWidth',1.2);         % boundary
hold on;% 'BinWidth',0.1, ...
histogram(y3, 15, 'Normalization','probability',...
              'EdgeColor','r', ...      % boundary
              'FaceColor','none', ...   % fill
              'LineWidth',1.2);         % boundary
ylim([0 0.6]);
yt = get(gca,'YTick');            % 
set(gca, 'YTick', yt, 'YTickLabel', yt*100);

xlabel('U_{BSE} (ppb)');ylabel('Frequency (%)');
title(['Histogram of Solutions Satisfying ',resid_col_names{1}], 'FontSize', 12,'FontWeight', 'bold');
set(gca,'FontUnits','points','FontSize',7,'FontWeight','bold',...
    'Linewidth',1.5,'LabelFontSizeMultiplier',1.0);
text(-0.126,1.01,'a','Units','normalized','FontWeight','bold','FontSize',8);

gca=subplot(3,2,2);
% histogram of Rs
k=4;
y1=success_F(:,k)/1e22;
y2=success_T(:,k)/1e22;
y3=success_Xe(:,k)/1e22;

histogram(y1, 15, 'Normalization','probability',...            
              'EdgeColor','k', ...      % boundary
              'FaceColor','none', ...   % fill
              'LineWidth',1.2);         % boundary
hold on;% 'BinWidth',0.1, ...
histogram(y2, 15, 'Normalization','probability',...
              'EdgeColor','b', ...      % boundary
              'FaceColor','none', ...   % fill
              'LineWidth',1.2);         % boundary
histogram(y3, 15, 'Normalization','probability',...
              'EdgeColor','r', ...      % boundary
              'FaceColor','none', ...   % fill
              'LineWidth',1.2);         % boundary

ylim([0 0.6]);
yt = get(gca,'YTick');            % 
set(gca, 'YTick', yt, 'YTickLabel', yt*100);

           % 
xlabel('R_{s} (e22 kg/Gyr)');ylabel('Frequency (%)');
set(gca,'FontUnits','points','FontSize',7,'FontWeight','bold',...
    'Linewidth',1.5,'LabelFontSizeMultiplier',1.0);
text(-0.126,1.01,'b','Units','normalized','FontWeight','bold','FontSize',8);

h=legend('CC','Thermal','Xe');
set(h,'Orientation','horizontal','Box','off',...
    'ItemTokenSize',[18 18],...
    'Position',[0.56 0.87 0.58 0.035],'FontSize',6);


% for 128Xe/130Xe
success_mask = abs( results(:, Erxe128v130_idx) ) <= 1;
success_Xe=results(success_mask,:);

gca=subplot(3,2,3);
% histogram of Ubse
k=6;
y2=success_T(:,k)*1e9;
y3=success_Xe(:,k)*1e9;


histogram(y2, 15, 'Normalization','probability',...            
              'EdgeColor','b', ...      % boundary
              'FaceColor','none', ...   % fill
              'LineWidth',1.2);         % boundary
hold on;% 'BinWidth',0.1, ...
histogram(y3, 15, 'Normalization','probability',...
              'EdgeColor','r', ...      % boundary
              'FaceColor','none', ...   % fill
              'LineWidth',1.2);         % boundary
ylim([0 0.6]);
yt = get(gca,'YTick');            % 
set(gca, 'YTick', yt, 'YTickLabel', yt*100);

xlabel('U_{BSE} (ppb)');ylabel('Frequency (%)');
title(['Histogram of Solutions Satisfying ',resid_col_names{2}], 'FontSize', 12,'FontWeight', 'bold');
set(gca,'FontUnits','points','FontSize',7,'FontWeight','bold',...
    'Linewidth',1.5,'LabelFontSizeMultiplier',1.0);
text(-0.126,1.01,'c','Units','normalized','FontWeight','bold','FontSize',8);

gca=subplot(3,2,4);
% histogram of Rs
k=4;
y1=success_F(:,k)/1e22;
y2=success_T(:,k)/1e22;
y3=success_Xe(:,k)/1e22;

histogram(y1, 15, 'Normalization','probability',...            
              'EdgeColor','k', ...      % boundary
              'FaceColor','none', ...   % fill
              'LineWidth',1.2);         % boundary
hold on;% 'BinWidth',0.1, ...
histogram(y2, 15, 'Normalization','probability',...
              'EdgeColor','b', ...      % boundary
              'FaceColor','none', ...   % fill
              'LineWidth',1.2);         % boundary
histogram(y3, 15, 'Normalization','probability',...
              'EdgeColor','r', ...      % boundary
              'FaceColor','none', ...   % fill
              'LineWidth',1.2);         % boundary

ylim([0 0.6]);
yt = get(gca,'YTick');            % 
set(gca, 'YTick', yt, 'YTickLabel', yt*100);
text(-0.126,1.01,'d','Units','normalized','FontWeight','bold','FontSize',8);

           % 
xlabel('R_{s} (e22 kg/Gyr)');ylabel('Frequency (%)');
set(gca,'FontUnits','points','FontSize',7,'FontWeight','bold',...
    'Linewidth',1.5,'LabelFontSizeMultiplier',1.0);


% for 136Xe/132Xe
success_mask = abs( results(:, Erxe136v132_idx) ) <= 1;
success_Xe=results(success_mask,:);
gca=subplot(3,2,5);
% histogram of Ubse
k=6;
y2=success_T(:,k)*1e9;
y3=success_Xe(:,k)*1e9;


histogram(y2, 15, 'Normalization','probability',...            
              'EdgeColor','b', ...      % boundary
              'FaceColor','none', ...   % fill
              'LineWidth',1.2);         % boundary
hold on;% 'BinWidth',0.1, ...
histogram(y3, 15, 'Normalization','probability',...
              'EdgeColor','r', ...      % boundary
              'FaceColor','none', ...   % fill
              'LineWidth',1.2);         % boundary
ylim([0 0.6]);
yt = get(gca,'YTick');            % 
set(gca, 'YTick', yt, 'YTickLabel', yt*100);

xlabel('U_{BSE} (ppb)');ylabel('Frequency (%)');
title(['Histogram of Solutions Satisfying ',resid_col_names{6}], 'FontSize', 12,'FontWeight', 'bold');
set(gca,'FontUnits','points','FontSize',7,'FontWeight','bold',...
    'Linewidth',1.5,'LabelFontSizeMultiplier',1.0);
text(-0.126,1.01,'e','Units','normalized','FontWeight','bold','FontSize',8);

gca=subplot(3,2,6);
% histogram of Rs
k=4;
y1=success_F(:,k)/1e22;
y2=success_T(:,k)/1e22;
y3=success_Xe(:,k)/1e22;

histogram(y1, 15, 'Normalization','probability',...            
              'EdgeColor','k', ...      % boundary
              'FaceColor','none', ...   % fill
              'LineWidth',1.2);         % boundary
hold on;% 'BinWidth',0.1, ...
histogram(y2, 15, 'Normalization','probability',...
              'EdgeColor','b', ...      % boundary
              'FaceColor','none', ...   % fill
              'LineWidth',1.2);         % boundary
histogram(y3, 15, 'Normalization','probability',...
              'EdgeColor','r', ...      % boundary
              'FaceColor','none', ...   % fill
              'LineWidth',1.2);         % boundary

ylim([0 0.6]);
yt = get(gca,'YTick');            % 
set(gca, 'YTick', yt, 'YTickLabel', yt*100);

           % 
xlabel('R_{s} (e22 kg/Gyr)');ylabel('Frequency (%)');
set(gca,'FontUnits','points','FontSize',7,'FontWeight','bold',...
    'Linewidth',1.5,'LabelFontSizeMultiplier',1.0);
text(-0.126,1.01,'f','Units','normalized','FontWeight','bold','FontSize',8);



%% % S4 histogram: Xei=3.2e8

% load data
success_F_struct=load('cc/out_cc_grwoth_success_lhs_MC_Krw_lowf.mat');
success_F=success_F_struct.success_F;
[size_F,nparam_cc]=size(success_F);

success_T_struct=load('thermal/out_thermal_success_lhs_test_fixK_U_MC_lowerT_Rac1100_lowf.mat');
success_T=success_T_struct.success_T;
[size_T,nparam_T]=size(success_T);

results_Xe_struct=load('Xe/out_Xe_lhs_new_3e8_test_fixK_U_MC_lowerT_alot_Rac1100_new_lowf_cc1.mat');
results=results_Xe_struct.results;

col_names = {'ts_cc','kappa_gcc','kappa_rcc','Rs','Rp',...
    'Ubse','Qc_pd', 'dQc', 'Ti', 'eta_ref',...
    'Xm_init','Ms_init','Fr_w',...
    'Xe_init','Fr_Xe',...
    'Fd_mor','Fd_p',...
    'Er_Pmor_pd',... % 'RMSE_mor_processing rate',...
    'Er_xec','Er_xe128v130','Er_xe128v132','Er_xe130v132','Er_xe131v132','Er_xe134v132','Er_xe136v132'};

param_names = {'ts','\kappa_g','\kappa_r','Rs','Rp',...
    'Ubse','Qc_pd', 'Ti', 'eta_ref',...
    'Fr_Xe',...
    'Fd_mor','Fd_p'};

Erxec_idx = find(strcmp(col_names, 'Er_xec'));
Erxe128v130_idx = find(strcmp(col_names, 'Er_xe128v130'));
Erxe128v132_idx = find(strcmp(col_names, 'Er_xe128v132'));
Erxe130v132_idx = find(strcmp(col_names, 'Er_xe130v132'));
Erxe131v132_idx = find(strcmp(col_names, 'Er_xe131v132'));
Erxe134v132_idx = find(strcmp(col_names, 'Er_xe134v132'));
Erxe136v132_idx = find(strcmp(col_names, 'Er_xe136v132'));

resid_col_names = {'^{130}Xe','^{128}Xe/^{130}Xe','^{130}Xe/^{132}Xe','^{131}Xe/^{132}Xe','^{134}Xe/^{132}Xe','^{136}Xe/^{132}Xe'};
resid_col_ids = [Erxec_idx,Erxe128v130_idx, Erxe130v132_idx, Erxe131v132_idx, Erxe134v132_idx, Erxe136v132_idx];


figure;
set(gcf,'color','white','Units','centimeters','Position',[2 5 17.8 12]);

% for Xe concentration
success_mask = abs( results(:, Erxec_idx) ) <= 1;
success_Xe=results(success_mask,:);
gca=subplot(3,2,1);
% histogram of Ubse
k=6;
y2=success_T(:,k)*1e9;
y3=success_Xe(:,k)*1e9;


histogram(y2, 15, 'Normalization','probability',...            
              'EdgeColor','b', ...      % boundary
              'FaceColor','none', ...   % fill
              'LineWidth',1.2);         % boundary
hold on;% 'BinWidth',0.1, ...
histogram(y3, 15, 'Normalization','probability',...
              'EdgeColor','r', ...      % boundary
              'FaceColor','none', ...   % fill
              'LineWidth',1.2);         % boundary
ylim([0 0.6]);
yt = get(gca,'YTick');            % 
set(gca, 'YTick', yt, 'YTickLabel', yt*100);

xlabel('U_{BSE} (ppb)');ylabel('Frequency (%)');
title(['Histogram of Solutions Satisfying ',resid_col_names{1}], 'FontSize', 12,'FontWeight', 'bold');
set(gca,'FontUnits','points','FontSize',7,'FontWeight','bold',...
    'Linewidth',1.5,'LabelFontSizeMultiplier',1.0);
text(-0.126,1.01,'a','Units','normalized','FontWeight','bold','FontSize',8);

gca=subplot(3,2,2);
% histogram of Rs
k=4;
y1=success_F(:,k)/1e22;
y2=success_T(:,k)/1e22;
y3=success_Xe(:,k)/1e22;

histogram(y1, 15, 'Normalization','probability',...            
              'EdgeColor','k', ...      % boundary
              'FaceColor','none', ...   % fill
              'LineWidth',1.2);         % boundary
hold on;% 'BinWidth',0.1, ...
histogram(y2, 15, 'Normalization','probability',...
              'EdgeColor','b', ...      % boundary
              'FaceColor','none', ...   % fill
              'LineWidth',1.2);         % boundary
histogram(y3, 15, 'Normalization','probability',...
              'EdgeColor','r', ...      % boundary
              'FaceColor','none', ...   % fill
              'LineWidth',1.2);         % boundary

ylim([0 0.6]);
yt = get(gca,'YTick');            % 
set(gca, 'YTick', yt, 'YTickLabel', yt*100);

           % 
xlabel('R_{s} (e22 kg/Gyr)');ylabel('Frequency (%)');
set(gca,'FontUnits','points','FontSize',7,'FontWeight','bold',...
    'Linewidth',1.5,'LabelFontSizeMultiplier',1.0);
text(-0.126,1.01,'b','Units','normalized','FontWeight','bold','FontSize',8);

h=legend('CC','Thermal','Xe');
set(h,'Orientation','horizontal','Box','off',...
    'ItemTokenSize',[18 18],...
    'Position',[0.56 0.87 0.58 0.035],'FontSize',6);

% for 128Xe/130Xe
success_mask = abs( results(:, Erxe128v130_idx) ) <= 1;
success_Xe=results(success_mask,:);

gca=subplot(3,2,3);
% histogram of Ubse
k=6;
y2=success_T(:,k)*1e9;
y3=success_Xe(:,k)*1e9;


histogram(y2, 15, 'Normalization','probability',...            
              'EdgeColor','b', ...      % boundary
              'FaceColor','none', ...   % fill
              'LineWidth',1.2);         % boundary
hold on;% 'BinWidth',0.1, ...
histogram(y3, 15, 'Normalization','probability',...
              'EdgeColor','r', ...      % boundary
              'FaceColor','none', ...   % fill
              'LineWidth',1.2);         % boundary
ylim([0 0.6]);
yt = get(gca,'YTick');            % 
set(gca, 'YTick', yt, 'YTickLabel', yt*100);

xlabel('U_{BSE} (ppb)');ylabel('Frequency (%)');
title(['Histogram of Solutions Satisfying ',resid_col_names{2}], 'FontSize', 12,'FontWeight', 'bold');
set(gca,'FontUnits','points','FontSize',7,'FontWeight','bold',...
    'Linewidth',1.5,'LabelFontSizeMultiplier',1.0);
text(-0.126,1.01,'c','Units','normalized','FontWeight','bold','FontSize',8);

gca=subplot(3,2,4);
% histogram of Rs
k=4;
y1=success_F(:,k)/1e22;
y2=success_T(:,k)/1e22;
y3=success_Xe(:,k)/1e22;

histogram(y1, 15, 'Normalization','probability',...            
              'EdgeColor','k', ...      % boundary
              'FaceColor','none', ...   % fill
              'LineWidth',1.2);         % boundary
hold on;% 'BinWidth',0.1, ...
histogram(y2, 15, 'Normalization','probability',...
              'EdgeColor','b', ...      % boundary
              'FaceColor','none', ...   % fill
              'LineWidth',1.2);         % boundary
histogram(y3, 15, 'Normalization','probability',...
              'EdgeColor','r', ...      % boundary
              'FaceColor','none', ...   % fill
              'LineWidth',1.2);         % boundary

ylim([0 0.6]);
yt = get(gca,'YTick');            % 
set(gca, 'YTick', yt, 'YTickLabel', yt*100);
text(-0.126,1.01,'d','Units','normalized','FontWeight','bold','FontSize',8);

           % 
xlabel('R_{s} (e22 kg/Gyr)');ylabel('Frequency (%)');
set(gca,'FontUnits','points','FontSize',7,'FontWeight','bold',...
    'Linewidth',1.5,'LabelFontSizeMultiplier',1.0);


% for 136Xe/132Xe
success_mask = abs( results(:, Erxe136v132_idx) ) <= 1;
success_Xe=results(success_mask,:);
gca=subplot(3,2,5);
% histogram of Ubse
k=6;
y2=success_T(:,k)*1e9;
y3=success_Xe(:,k)*1e9;


histogram(y2, 15, 'Normalization','probability',...            
              'EdgeColor','b', ...      % boundary
              'FaceColor','none', ...   % fill
              'LineWidth',1.2);         % boundary
hold on;% 'BinWidth',0.1, ...
histogram(y3, 15, 'Normalization','probability',...
              'EdgeColor','r', ...      % boundary
              'FaceColor','none', ...   % fill
              'LineWidth',1.2);         % boundary
ylim([0 0.6]);
yt = get(gca,'YTick');            % 
set(gca, 'YTick', yt, 'YTickLabel', yt*100);

xlabel('U_{BSE} (ppb)');ylabel('Frequency (%)');
title(['Histogram of Solutions Satisfying ',resid_col_names{6}], 'FontSize', 12,'FontWeight', 'bold');
set(gca,'FontUnits','points','FontSize',7,'FontWeight','bold',...
    'Linewidth',1.5,'LabelFontSizeMultiplier',1.0);
text(-0.126,1.01,'e','Units','normalized','FontWeight','bold','FontSize',8);

gca=subplot(3,2,6);
% histogram of Rs
k=4;
y1=success_F(:,k)/1e22;
y2=success_T(:,k)/1e22;
y3=success_Xe(:,k)/1e22;

histogram(y1, 15, 'Normalization','probability',...            
              'EdgeColor','k', ...      % boundary
              'FaceColor','none', ...   % fill
              'LineWidth',1.2);         % boundary
hold on;% 'BinWidth',0.1, ...
histogram(y2, 15, 'Normalization','probability',...
              'EdgeColor','b', ...      % boundary
              'FaceColor','none', ...   % fill
              'LineWidth',1.2);         % boundary
histogram(y3, 15, 'Normalization','probability',...
              'EdgeColor','r', ...      % boundary
              'FaceColor','none', ...   % fill
              'LineWidth',1.2);         % boundary

ylim([0 0.6]);
yt = get(gca,'YTick');            % 
set(gca, 'YTick', yt, 'YTickLabel', yt*100);

           % 
xlabel('R_{s} (e22 kg/Gyr)');ylabel('Frequency (%)');
set(gca,'FontUnits','points','FontSize',7,'FontWeight','bold',...
    'Linewidth',1.5,'LabelFontSizeMultiplier',1.0);
text(-0.126,1.01,'f','Units','normalized','FontWeight','bold','FontSize',8);


exportgraphics(gcf, 'test_figS4_histogram_3e8.pdf', ...
    'ContentType','image', ...   % vector/image
    'BackgroundColor','white', ... % 
    'Resolution',600); 


%% % S5 shap value

% load data
data = load('Xe/SHAP_interactions.mat');
SHAP_inter = data.SHAP_interactions;  % (n_targets, n_features, n_features)
feature_names = {'t_{s}','\kappa_g','\kappa_r','R_s','R_p',...
    'U_{BSE}','Q_c', 'T_i', '\eta_{ref}',...
    'Fr',...
    'Fd_{M}','Fd_{P}'};
target_names = {'^{130}Xe','^{128}Xe/^{130}Xe',...
    '^{130}Xe/^{132}Xe','^{131}Xe/^{132}Xe',...
    '^{134}Xe/^{132}Xe','^{136}Xe/^{132}Xe'};

n_targets = size(SHAP_inter, 1);
n_features = size(SHAP_inter, 2);

% 
top_k = 5;  % 

figure;
set(gcf,'color','white','Units','centimeters','Position',[2 5 17.8 12]);

for target_idx = 1:n_targets
    % interactions for each target
    inter_matrix = squeeze(SHAP_inter(target_idx, :, :));
    
    % only upper to avoid repeating
    upper_tri_indices = triu(true(n_features), 1);
    inter_values = inter_matrix(upper_tri_indices);
    
    % get index
    [row_idx, col_idx] = find(upper_tri_indices);
    
    % sort
    [sorted_values, sort_idx] = sort(inter_values, 'descend');
    top_indices = sort_idx(1:min(top_k, length(sort_idx)));
    
    % 
    top_values = sorted_values(1:min(top_k, length(sort_idx)));
    top_row_idx = row_idx(top_indices);
    top_col_idx = col_idx(top_indices);
    
    % mark the name for each interaction
    labels = cell(length(top_values), 1);
    for k = 1:length(top_values)
        feat1 = feature_names{top_row_idx(k)};
        feat2 = feature_names{top_col_idx(k)};
        labels{k} = sprintf('%s x %s', feat1, feat2);
    end
    
    % print
    fprintf('\n%s - Top %d interactions:\n', target_names{target_idx}, top_k);
    for k = 1:length(top_values)
        fprintf('  %d. %s: %.4e\n', k, labels{k}, top_values(k));
    end
    
    % ====================================================================
    % 3. plot
    % ====================================================================
    gca=subplot(ceil(n_targets/2), 2, target_idx);
    
    %  largest interactions at the top
    flipped_values = flip(top_values);
    flipped_labels = flip(labels);
    
    % 
    barh(flipped_values, 'FaceColor', [0.26, 0.45, 0.64]);
    
    % 
    set(gca, 'YTick', 1:length(flipped_labels));
    set(gca, 'YTickLabel', flipped_labels);
    
    % 
    set(gca, 'YTick', 1:length(labels));
    set(gca, 'YTickLabel', flip(labels));
    
    % 
    title(target_names{target_idx}, 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter','tex');
    % xlabel('mean |SHAP interaction|', 'FontSize', 11);
    
    % 
    grid on;
    set(gca, 'FontSize', 10);
    xlim([0, max(top_values)*1.1]);  % 
    grid on;
    set(gca,'FontUnits','points','FontSize',7,'FontWeight','bold',...
    'Linewidth',1.5,'LabelFontSizeMultiplier',1.0);
    % set(gca,'Position',[0.06 0.7177 0.505 0.2]);
end

sgtitle('Top 5 SHAP Feature Interactions', 'FontSize', 16, 'FontWeight', 'bold');

% exportgraphics(gcf, 'test_figS5_shap.pdf', ...
%     'ContentType','image', ...   % vector/image
%     'BackgroundColor','white', ... % 
%     'Resolution',600); 


