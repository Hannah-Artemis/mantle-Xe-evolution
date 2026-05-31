%% Random Forest on Xe - fin

% load data
results_Xe_struct=load('Xe/out_Xe_lhs_new_3e7_test_fixK_U_MC_lowerT_alot.mat');
results=results_Xe_struct.results;
[size_Xe,nparam_Xe]=size(results);

% Get parameters and RMSE
X_rf = results(:, [1:7,9:10,15:17]);  % 12 params

col_names = {'ts_cc','kappa_gcc','kappa_rcc','Rs','Rp',...
    'Ubse','Qc_pd', 'dQc', 'Ti', 'eta_ref',...
    'Xm_init','Ms_init','Fr_w',...
    'Xe_init','Fr_Xe',...
    'Fd_mor','Fd_p',...
    'Er_Pmor_pd',... % 'RMSE_water',...
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


% 
% for 3e6 samples, each model needs 700s
% 6 models --> 4800s
RF_models = cell(length(resid_col_ids), 1);
RF_importance = zeros(length(resid_col_ids), length(param_names));

tic;
for i = 1:length(resid_col_ids)

    disp(resid_col_names{i});
    
    resid_col = resid_col_ids(i);
    rmse = abs(results(:, resid_col));  % 
    % error = results(:, resid_col);

    % Define success based on RMSE threshold
    threshold = 1;
    Y = double( rmse < threshold);  % success = 1, failure = 0


    % train RF
    mdl = fitrensemble(X_rf, Y, ...
        'Method','Bag', ...
        'NumLearningCycles',100, ...
        'Learners','tree');
    
    % 
    RF_models{i} = mdl;    
    imp = predictorImportance(mdl);
    imp_percent = 100 * imp / sum(imp);
    RF_importance(i, :) = imp_percent;


    % plot
    figure;
    set(gcf,'color','white');
    bar(RF_importance(i, :)); 
    xticks(1:length(param_names));
    set(gca, 'XTickLabel', param_names);
    ylabel('Feature Importance');
    title(['Influence on ', resid_col_names{i}], 'FontSize', 16,'FontWeight', 'bold');
    grid on;

    elapsedTime = toc;  % s
    disp(['used time = ', num2str(elapsedTime), ' s']);

end


% % save
% save('Xe/Xe_rf_residuals_Xei3p2e7_alot.mat', 'RF_models', 'RF_importance', ...
%      'resid_col_names', 'param_names','-v7.3');


