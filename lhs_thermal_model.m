%%

% LHS for thermal model
% based on successful [kappa_gcc,ts] model from stage1- cc growth

clear all;
tic;
%
ncores=6;
if isempty(gcp('nocreate'))
    parpool('local', ncores);  % 
end

par_solver_new;

N = 1e5; % Number of samples
Nseeds_stage1 = 50;
%Nseeds_stage1 = 1;

% Open output file
%fileID = fopen('out_phase1a_lhs.dat','w');

%% Set the time period
tmax = t_pd/yr_s/1e9;% the age of solar system, in unit Ga
dt = 0.001;% length of each timestep
t = 0:dt:tmax;
nt = length(t);% number of timesteps
t = t';

% default_input_for_degassing;
% t=Simulation_input.time_series;
% t=t/yr_s/1e9;
% nt = length(t);

%% Define the a priori ranges for independent variables
%% Step 1: Latin Hypercube Sampling (LHS)
%% Constants used in the model
Mcp = mcc_val;% mass of continental crust at present-day, in unit kg
eta_ref = 7.5e14;% fixed eta_ref from the relationship of present day Tp & Qs
para_Scale=2;% 2: Schubert
Hsf=2;

%load successful cc growth from stage1  
success_cc_struct=load('cc/out_cc_grwoth_success_lhs_MC_Krw.mat');
success_cc=success_cc_struct.success_F;
% seeds_stage2_struct=load('cc/out_for_seeds_stage2.mat');
% seeds_stage2=seeds_stage2_struct.seeds_stage2;
% 

[Ncc,nparam]=size(success_cc);


%col_names = {'ts','kappa_g', 'kappa_r', 'Rs', 'Rp', 'RMSE'};

% Define parameter bounds

param_bounds = [
   10e-9, 20e-9;     % BSE type or Ubse
%     5e12, 15e12;    % Qcpd: present day core heat flux [W]
%     2e12, 5e12;   % dQc: variation of core heat flux
    5e12, 15e12;    % Qcpd: present day core heat flux [W]
    0, 0;   % dQc: variation of core heat flux
    2200, 2700;   % Ti: initial mantle average temperature [2200 2800] [2300 2700]
    % 5e14, 19e14    % eta_ref: constrained by the relationship of today Qs & Tp                 
    7.39e14, 2.5e15    % eta_ref: constrained by the relationship of today Qs & Tp                 
];

% LHS sampling for continous value
n_params = 5; % number of model parameters BSE, Qcpd, dQc, Ti)
lhs_unit = lhsdesign(N, n_params);
params = zeros(N, n_params);
for i = 1:n_params
    params(:, i) = param_bounds(i,1) + (param_bounds(i,2)-param_bounds(i,1)) * lhs_unit(:, i);
end

% % BSE
% discrete_values = [-5, -4, -3, -2, -1, 1, 2, 3, 4, 5];
% % discrete_values = [-4, -1, 1, 2, 3, 4, 5]; % only fixed K/U 
% % discrete_values = [ 2 ]; % only fixed K/U 
% %  discrete_values
% idxs = ceil(lhs_unit(:,1) * numel(discrete_values)); % random number [1,6]
% idxs(idxs < 1) = 1;  % lower boundary
% idxs(idxs > numel(discrete_values)) = numel(discrete_values); %upper boundary
% params(:,1) = discrete_values(idxs);


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

%% Submit parfeval jobs
% futures(N*Nseeds_stage1) = parallel.FevalFuture;
futures(N) = parallel.FevalFuture;
k=0;
for i = 1:N
    index_cc=randi(Ncc);
    j=index_cc;
   % for j = 1:Nseeds_stage1
        k=k+1;
        disp(k);
        param_struct = struct('time_cc_input', success_cc(j,1)*1e9*yr_s,...% seeds_stage2(j,1)*1e9*yr_s, ... %change from Ga to s
                              'kappa_gcc', success_cc(j,2),...% seeds_stage2(j,2), ...
                              'kappa_rcc', success_cc(j,3),...% seeds_stage2(j,3), ...
                              'Rs', success_cc(j,4),...% seeds_stage2(j,4), ...
                              'Rp', success_cc(j,5),...% seeds_stage2(j,5), ...
                              'BSE',      -4, ...
                              'Ubse', params(i,1),...
                              'Qc_pd_input',      params(i,2), ...
                              'dQc_input',      params(i,3), ...
                              'Ti',        params(i,4),...
                              'Mcp',     Mcp, ...
                              'eta_ref',        params(i,5), ...
                              'para_Scale',    para_Scale,  ...
                              'Hsf',      Hsf, ...
                              't',       t*1e9*yr_s, ... %change from Ga to s
                              'tmax',    tmax*1e9*yr_s, ...
                              'dt',     dt,...  % in Ga, used only for get index for Tp_harz at given time in Ga
                               'Tp_obs1', Tp_obs1, ...
                              'Tp_obs2', Tp_obs2,...
                              'Qs_obs', Qs_range);
        futures(k) = parfeval(@run_thermal_model, 1, param_struct);
  %  end
end


%% Collect outputs
% results = zeros(N*Nseeds_stage1,13);
results = zeros(N,13);
for i = 1:N %*Nseeds_stage1
    try
        results(i, :) = fetchOutputs(futures(i));
    catch ME
        warning('Error in task %d: %s', i, ME.message);
        results(i, :) = nan(1,13);  % 
    end
end

%% Save results
results(:,1)=results(:,1)/1e9/yr_s;
col_names = {'ts_cc','kappa_gcc','kappa_rcc','Rs','Rp','Ubse','Qc_pd', 'dQc', 'Ti', 'eta_ref','RMSE_pd','RMSE_t','RMSE_Q'};
% save('thermal/out_thermal_lhs_test_MC_lowerT_Rac1100.mat', 'results', 'col_names','-v7.3');

% delete(gcp('nocreate'));

elapsedTime = toc;  % s
disp(['used time = ', num2str(elapsedTime), ' s']);



%% postprocesses
%% select successful solutions for plotting results
crit_misfit_T1 = 0.5; %1:[1250 1450];0.5:[1300 1400];0.2:[1330 1370];
crit_misfit_T2 = 1;
crit_misfit_Q = 1; % within the Qs range

rmse_idx = find(strcmp(col_names, 'RMSE_pd'));
rmset_idx = find(strcmp(col_names, 'RMSE_t'));
rmseq_idx = find(strcmp(col_names, 'RMSE_Q'));
success_mask = ( results(:, rmse_idx) < crit_misfit_T1 & results(:, rmset_idx) < crit_misfit_T2 & results(:, rmseq_idx) < crit_misfit_Q );
success_T = results(success_mask, :);

% whether sampling number is enough
if size(success_T, 1) < 1000
    disp('Success sample not enough');
end
% save('thermal/out_thermal_success_lhs_test_fixK_U_MC_lowerT_Rac1100.mat', 'success_T', '-v7.3');


%% Calculate Temperature curve of successful models
% success_T_struct=load('thermal/out_thermal_success_lhs_test_fixK_U.mat');
% success_T=success_T_struct.success_T;

[size_T,nparameter]= size(success_T);
NGcc = nan(nt,size_T);
Rcc = nan(nt,size_T);
Gcc = nan(nt,size_T);
Tm = nan(nt,size_T);
Tp = nan(nt,size_T);
Qs = nan(nt,size_T);
Qc = nan(nt,size_T);
H = nan(nt,size_T);
Ra = nan(nt,size_T);
for l=1:size_T
    if mod(l,10)==0
        disp(['l=' num2str(l) ' of ' num2str(size_T)]);% keep track of the calculation
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
    [NGcc_model,Rcc_model,Gcc_model] = CC_growth_fun1(t,ts_model,tmax,...
    Mcp,kappa_g_model,Rp_model,Rs_model,kappa_r_model);
    % save the results
    NGcc(:,l) = NGcc_model;
    Rcc(:,l) = Rcc_model;
    Gcc(:,l) = Gcc_model;

    % Calculate the corresponding Temperature curve
    [Tp_model,Tm_model,Qs_model,Qc_model,H_model,Ra_model,eta_out] = Thermal_model(t*1e9*yr_s,Ti_model,eta_model,Qc_model,dQc_model,...
            BSE_model,Ubse_model,kappa_g_model,ts_model*1e9*yr_s,para_Scale,Hsf);
    Tp(:,l) = Tp_model;
    Tm(:,l) = Tm_model;
    Qs(:,l) = Qs_model;
    Qc(:,l) = Qc_model;
    H(:,l) = H_model;
    Ra(:,l) = Ra_model;

end % for l=1:itermax


%%
[NGcc_5,NGcc_25,NGcc_50,NGcc_75,NGcc_95] = calculate_percentile_fun(NGcc,nt,t,size_T);
[Rcc_5,Rcc_25,Rcc_50,Rcc_75,Rcc_95] = calculate_percentile_fun(Rcc,nt,t,size_T);
[Gcc_5,Gcc_25,Gcc_50,Gcc_75,Gcc_95] = calculate_percentile_fun(Gcc,nt,t,size_T);
[Tp_5,Tp_25,Tp_50,Tp_75,Tp_95] = calculate_percentile_fun(Tp,nt,t,size_T);
[Qs_5,Qs_25,Qs_50,Qs_75,Qs_95] = calculate_percentile_fun(Qs/1e12,nt,t,size_T);

figure(4);
subplot(2,2,1); hold off;
plot(t,NGcc_25,'k--',t,NGcc_75,'k--',t,NGcc_5,'k:',t,NGcc_95,'k:',t,NGcc_50,'r-');
xlabel('Time (Gyr)'); ylabel('Net crustal growth (kg)'); grid on;
subplot(2,2,2); hold off;
hold off;
semilogy(t,Gcc_25,'k--',t,Gcc_75,'k--',t,Gcc_5,'k:',t,Gcc_95,'k:',t,Gcc_50,'r-'); hold on;
xlabel('Time (Gyr)'); ylabel('Crustal generation rate (kg/Ga)'); grid on;
subplot(2,2,3); hold off;
plot(t,Qs_25,'k--',t,Qs_75,'k--',t,Qs_5,'k:',t,Qs_95,'k:',t,Qs_50,'r-'); hold on;
xlabel('Time (Gyr)'); ylabel('Surface Heat Flow (TW)'); grid on;
subplot(2,2,4); hold off;
plot(t,Tp_25,'k--',t,Tp_75,'k--',t,Tp_5,'k:',t,Tp_95,'k:',t,Tp_50,'r-');
xlabel('Time (Gyr)'); ylabel('Potential Temperature (°C)'); grid on;

[Tp_lower,Tp_upper] = calculate_margin_fun(Tp,nt,t,size_T);
[Qs_lower,Qs_upper] = calculate_margin_fun(Qs,nt,t,size_T);

color_unit=[0.6, 0.8, 1.0];
color_purple = [0.5, 0.2, 0.7];

figure(502);
t = t(:)';                        
Tp_lower = Tp_lower(:)';      
Tp_upper = Tp_upper(:)';
Tp_5 = Tp_5(:)';      
Tp_95 = Tp_95(:)';
Tp_25 = Tp_25(:)';      
Tp_75 = Tp_75(:)';
set(gcf,'color','white');
x_fill = [t, fliplr(t)];
y_fill1 = [Tp_5,(fliplr(Tp_95))];
y_fill2 = [Tp_25,(fliplr(Tp_75))];
y_fill0 = [Tp_lower,(fliplr(Tp_upper))];
fill(x_fill, y_fill1,color_purple, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on;
fill(x_fill, y_fill2,color_purple, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on;
%fill(x_fill, y_fill0,color_purple, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on;

plot(t, Tp_50, 'm-', 'LineWidth', 1.5);
scatter( t_Herz,Tp_Herz, 100, ...
    'o', 'MarkerEdgeColor', [1, 0.5, 0], ...
    'MarkerFaceColor', [1, 0.3, 0], 'LineWidth', 1.5);
% % 
% stairs(t, F_50, 'Color', [1 0.4 0.2], 'LineWidth', 2); 
%
xlabel('Time (Ga)');
ylabel('Mantle Potential Temperature (°C)');
grid on;
xlim([0 4.6]);
%ylim([0,1]);
legend('accepted models','','median model','observations','location','best');
%text(2.1, 0.6, 'Korenaga (2018)', 'FontSize', 10, 'Rotation', 0);
%annotation('textarrow', [0.55 0.48], [0.55 0.65], 'String', '');


%%

% 


% %% clustering
% 
% 
% % % Ubse param instead of BSE types
% n_seeds =100;
% seeds_stage3 = zeros(n_seeds, size(success_T,2)); 
% 
% X_T = [success_T(:, 1:7),success_T(:, 9:10)];
% X_T_norm = normalize(X_T, 'zscore');  %
% mu = mean(X_T);
% sigma = std(X_T);
% [idx, C] = kmeans(X_T_norm, n_seeds); % class & center
% C_original = C .* sigma + mu;
% 
% for i = 1:n_seeds
%     class_i = find(idx == i);
%     if isempty(class_i), continue; end
%     dists = vecnorm(X_T_norm(class_i,:) - C(i,:), 2, 2);
%     [~, minidx] = min(dists); % each class, find the model closest to its center 
%     global_idx = class_i(minidx);
%     seeds_stage3(i, :) = success_T(global_idx, :);
% end
% 
% seeds_stage3_fin=seeds_stage3;
% % save('thermal/out_for_seeds_stage3_test_fixK_U_MC.mat', 'seeds_stage3_fin', '-v7.3');


% %% % discrete BSE type
% n_seed =10;  % 
% n_BSE=10;
% BSE_range=[-5, -4, -3, -2, -1, 1, 2, 3, 4, 5];
% 
% seeds_stage3 = zeros(n_seed*n_BSE, size(success_T,2));  
% % including other parameters & RMSE
% 
% k=0;
% for i=1:n_BSE
%     BSE = BSE_range(i);
%     sub_idx = find(success_T(:,6)==BSE);  %
%     X_T = success_T(sub_idx, 1:5);        % only CC params
%     if (sum(success_T(:,6)==BSE)==0), continue; end
%     X_T_norm = normalize(X_T, 'zscore');  %
%     mu = mean(X_T);
%     sigma = std(X_T);
%     [idx, C] = kmeans(X_T_norm, n_seed); % class & center
%     C_original = C .* sigma + mu;
%     for j = 1:n_seed        
%         class_j = find(idx == j);
%         if isempty(class_j), continue; end
%         dists = vecnorm(X_T_norm(class_j,:) - C(j,:), 2, 2);
%         [~, minidx] = min(dists); % each class, find the model closest to its center 
%         global_idx = sub_idx(class_j(minidx)); 
%         k=k+1;
%         % disp([i,j,k]);
%         seeds_stage3(k, :) = success_T(global_idx, :);
%     end
% end
% 
% disp(k);


% seeds_stage3_fin=seeds_stage3(1:k, :);



% save('thermal/out_for_seeds_stage3_test.mat', 'seeds_stage3_fin', '-v7.3');

% % ranking
% rank_idx='number';
% rmse_idx=11;
% switch rank_idx
%     case 'number'
%         % based on number of points
%         counts = histcounts(idx, 1:(n_seed+1));  % the number of points in each class
%         [~, rank_order] = sort(counts, 'descend');  % 
%     case 'rmse'
%         % based on RMSE
%         [~, best_idx] = min(success_T(:, rmse_idx));   % find the best point
%         best_point = success_T(best_idx, 1:5);         % 
%         center_dists = vecnorm(C - best_point, 2, 2);  % find the class closest to best point
%         [~, rank_order] = sort(center_dists);          % 
% end
% 
% seeds_stage3_sorted = seeds_stage3_fin(rank_order, :);
% save('thermal/out_for_seeds_stage3_sorted_number_test.mat', 'seeds_stage3_sorted', '-v7.3');
% 






