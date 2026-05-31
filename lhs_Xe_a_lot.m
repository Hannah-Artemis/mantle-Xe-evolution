%%

% LHS for Xe model with large runs 
% which needs to be handled by clusters, or the storage may be run out
% based on successful [kappa_gcc,ts] model from stage1- cc growth

clear all;

%
ncores=4;
if isempty(gcp('nocreate'))
    parpool('local', ncores);  % 
end

tic;


%% Submit parfeval jobs
N_all=6e3;
N=N_all;

N_one_time=2e3;

par_solver_new;


%% Set the time period, only used for CC model, not for degassing
tmax = t_pd/yr_s/1e9;% the age of solar system, in unit Ga

default_input_for_degassing;

t=Simulation_input.time_series;
t=t/yr_s/1e9;
nt = length(t);



%% Constants used in the model
Mcp = mcc_val;% mass of continental crust at present-day, in unit kg


Msi=1.5;% 

SwXe = logical(1);
Xei=3.2e8; % atoms/g



%% load successful seeds from stage 1-2 
% seeds_stage4_struct=load('water/out_for_seeds_stage4.mat');
% seeds_stage4=seeds_stage4_struct.seeds_stage4_sorted;

success_T_struct=load('thermal/out_thermal_success_lhs_test_fixK_U_MC_lowerT_Rac1100_lowf.mat');
%success_T_struct=load('thermal/out_thermal_success_lhs_test_fixK_U_MC_lowerT.mat');
success_T=success_T_struct.success_T;
[Nthermal,nparam]=size(success_T);

% seeds_stage4_struct=load('thermal/out_for_seeds_stage3_test_fixK_U_MC.mat');
% seeds_stage4=seeds_stage4_struct.seeds_stage3_fin;

%% Define the a priori ranges for independent variables
%% Step 1: Latin Hypercube Sampling (LHS)
% Define parameter bounds
param_bounds = [
   0, 0.5;     % FRXe
   0.5, 1; % Fd_mor
   0.1, 1; % Fd_p
   ];
%     3e7, 3e7; %Xe_init atoms/g
% ];

% LHS sampling for continous value
n_params =3; % number of model parameters BSE, Qcpd, dQc, Ti)
lhs_unit = lhsdesign(N, n_params);
params = zeros(N, n_params);
for i = 1:n_params
    params(:, i) = param_bounds(i,1) + (param_bounds(i,2)-param_bounds(i,1)) * lhs_unit(:, i);
end


%% load observed present day mantle processing rate
P_mor_pd = [1e14 10e14];%kg/yr


%% Load in the observed water & Xe observation
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




cluster_end=floor(N_all/N_one_time);

futures(N_one_time) = parallel.FevalFuture;
results = zeros(N_all,25);

for cluster=1:cluster_end
    Ns=1+(cluster-1)*N_one_time;
    Ne=cluster*N_one_time;
    % Nl=Ne-Ns+1;   
    k=0;
    for i = Ns:Ne
        index_thermal=randi(Nthermal);
        j=index_thermal;
        % for j = 1:Nseeds_stage4
            k=k+1;
            disp(i);
           % disp(max(t));
            param_struct = struct( 't',       t*1e9*yr_s, ... % change from Ga to s
                                  'tmax',    tmax*1e9*yr_s, ... % change from Ga to s ====== above: time vector =======
                                  'time_cc_input', success_T(j,1)*1e9*yr_s,... %seeds_stage4(j,1)*1e9*yr_s,... % 0.2*1e9*yr_s, ... %   0.2*1e9*yr_s, ... % change from Ga to s
                                  'kappa_gcc', success_T(j,2),...%seeds_stage4(j,2), ... %  6,... % 
                                  'kappa_rcc', success_T(j,3),...%seeds_stage4(j,3), ... %  1,... %
                                  'Rs', success_T(j,4),...%seeds_stage4(j,4), ... %  10e22, ... % 
                                  'Rp', success_T(j,5),...%seeds_stage4(j,5), ... %  1.5e22, ... % % ====== above: CC model =======
                                  'BSE',   1,...  % seeds_stage4(j,6), ...%  4,...% 
                                  'Ubse', success_T(j,6),...%seeds_stage4(j,6), ...% 
                                  'Qc_pd_input',       success_T(j,7),...%seeds_stage4(j,7), ... %  15e12, ... %
                                  'dQc_input',      0,...%seeds_stage4(j,8), ... % 5e12,... % 
                                  'Ti',        success_T(j,9),...% seeds_stage4(j,9),...    % 2500,...%                          
                                  'eta_ref',        success_T(j,10),...%seeds_stage4(j,10), ... %  7.5e14,...% % ====== above: Thermal model ====== 
                                  'Fd_mor',    params(i,2),...  % 0.5,...%seeds_stage4(j,16),... %1, ... % 
                                  'Fd_p',     params(i,3),...% 0.5,...%%seeds_stage4(j,17),... % 1, ... %  % ====== above: MOR & Plume degassing ====== 
                                  'P_mor_pd', P_mor_pd,...
                                  'Fr_w',     0,...% seeds_stage4(j,13),...  % 0.25, ...       
                                  'Xm_init',  0,...% seeds_stage4(j,11),... %  0, ... % 
                                  'Ms_init',    1.5,...% seeds_stage4(j,12),...
                                  'water_obs1', water_obs1,... % ====== above: watre model ====== 
                                  'SwXe',  SwXe,...
                                  'Fr_Xe',      params(i,1),...                             
                                  'Xe_init',    Xei,...                             
                                  'Xe_obs', Xe_obs ... % ====== above: Xe model ======
                                   );
    
            futures(k) = parfeval(@run_degassing_model, 1, param_struct);
            %out_test=run_degassing_model(param_struct);
        % end
    end

    %% Collect outputs
    % results = zeros(N*Nseeds_stage4,25);
    
    ik=0;
    for i = Ns:Ne %
        ik=ik+1;
        try
            results(i, :) = fetchOutputs(futures(ik));
        catch ME
            warning('Error in task %d: %s', ik, ME.message);
            results(i, :) = nan(1,25);  % 
        end
    end


end




elapsedTime = toc;  % s
disp(['used time = ', num2str(elapsedTime), ' s']);

results(:,1)=results(:,1)/1e9/yr_s; % time_cc from s to Ga
col_names = {'ts_cc','kappa_gcc','kappa_rcc','Rs','Rp',...
    'Ubse','Qc_pd', 'dQc', 'Ti', 'eta_ref',...
    'Xm_init','Ms_init','Fr_w',...
    'Xe_init','Fr_Xe',...
    'Fd_mor','Fd_p',...
    'Er_Pmor_pd',... % 'RMSE_water',...
    'Er_xec','Er_xe128v130','Er_xe128v132','Er_xe130v132','Er_xe131v132','Er_xe134v132','Er_xe136v132'};
%    'RMSE_xec','RMSE_xe128v130','RMSE_xe128v132','RMSE_xe130v132','RMSE_xe131v132','RMSE_xe134v132','RMSE_xe136v132'};
save('Xe/out_Xe_lhs_new_3e7_test_fixK_U_MC_lowerT_alot_Rac1100_new_lowf_cc1.mat', 'results', 'col_names','-v7.3');

%% postprocesses
%% select successful solutions for plotting results
crit_misfit_Pmor = 1;
crit_misfit_w = 1;
crit_misfit_xec = 1;
crit_misfit_xe128v130 = 1; 
crit_misfit_xe128v132 = 1; 
crit_misfit_xe130v132 = 1; 
crit_misfit_xe131v132 = 1; 
crit_misfit_xe134v132 = 1; 
crit_misfit_xe136v132 = 1; 

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
                    abs( results(:, Erxe130v132_idx) ) <= crit_misfit_xe130v132 & ...
                    abs( results(:, Erxe131v132_idx) ) <= crit_misfit_xe131v132 & ...
                    abs( results(:, Erxe134v132_idx) ) <= crit_misfit_xe134v132 & ...
                    abs( results(:, Erxe136v132_idx) ) <= crit_misfit_xe136v132 & ...
                    abs( results(:, ErPmor_idx) ) <= crit_misfit_Pmor ); %& ...
                    %results(:, rmsew_idx) <= crit_misfit_w & ...


% rmsew_idx = find(strcmp(col_names, 'RMSE_water'));
% rmsexec_idx = find(strcmp(col_names, 'RMSE_xec'));
% rmsexe128v130_idx = find(strcmp(col_names, 'RMSE_xe128v130'));
% rmsexe128v132_idx = find(strcmp(col_names, 'RMSE_xe128v132'));
% rmsexe130v132_idx = find(strcmp(col_names, 'RMSE_xe130v132'));
% rmsexe131v132_idx = find(strcmp(col_names, 'RMSE_xe131v132'));
% rmsexe134v132_idx = find(strcmp(col_names, 'RMSE_xe134v132'));
% rmsexe136v132_idx = find(strcmp(col_names, 'RMSE_xe136v132'));
% 
% success_mask = ( results(:, rmsew_idx) <= crit_misfit_w & ...
%                     results(:, rmsexec_idx) <= crit_misfit_xec & ...
%                     results(:, rmsexe128v130_idx) <= crit_misfit_xe128v130 & ... % results(:, rmsexe128v132_idx) <= crit_misfit_xe128v132 & ...                    
%                     results(:, rmsexe130v132_idx) <= crit_misfit_xe130v132 & ...
%                     results(:, rmsexe131v132_idx) <= crit_misfit_xe131v132 & ...
%                     results(:, rmsexe134v132_idx) <= crit_misfit_xe134v132 & ...
%                     results(:, rmsexe136v132_idx) <= crit_misfit_xe136v132 );

success_Xe = results(success_mask, :);

% whether sampling number is enough
if size(success_Xe, 1) < 10
    disp('Success sample not enough');
end
save('Xe/out_Xe_success_lhs_new_3e7_test_fixK_U_MC_negativFD_Rac1100_new_lowf_cc1.mat', 'success_Xe', '-v7.3');
