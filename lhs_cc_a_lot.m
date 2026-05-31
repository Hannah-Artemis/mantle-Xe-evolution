%%

% LHS for Xe model with large runs 
% which needs to be handled by clusters, or the storage may be run out
% based on successful [kappa_gcc,ts] model from stage1- cc growth

clear all;

%
ncores=40;
if isempty(gcp('nocreate'))
    parpool('local', ncores);  % 
end

tic;


%% Submit parfeval jobs
N_all=5e5;
N=N_all;

N_one_time=1e5;

par_solver_new;


n_params = 6; % number of model parameters (e.g., kappa_r, kappa_g, Rs, Rp, ts, Krw)

% Open output file
%fileID = fopen('out_phase1a_lhs.dat','w');

%% Set the time period
tmax = t_pd/yr_s/1e9;% the age of solar system, in unit Ga
dt = 0.001;% length of each timestep
t = 0:dt:tmax;
nt = length(t);% number of timesteps
t = t';


%% Define the a priori ranges for independent variables
%% Step 1: Latin Hypercube Sampling (LHS)
%% Constants used in the model
Mcp = mcc_val;% mass of continental crust at present-day, in unit kg

% Define parameter bounds
param_bounds = [
    -3, 3;     % kappa_r: decay constant for crustal recycling rate, in unit Gyr-1
    -1, 30;    % kappa_g: decay constant for crustal growth rate, in unit Gyr-1
    0, 1e23;   % Rs: the initial crustal recycling rate (t=ts), in unit kg/Gyr
    0, 2e22;   % Rp: the present-day crustal recycling rate (t=tp), in unit kg/Gyr
    t_pd/yr_s/1e9-4.51, t_pd/yr_s/1e9-4; % ts: the onset time for crustal growth and recycle, in unit Ga
    0.1, 0.8;  % frw: crustal reworking rate factor
];

% LHS sampling
lhs_unit = lhsdesign(N, n_params);
params = zeros(N, n_params);
for i = 1:n_params
    params(:, i) = param_bounds(i,1) + (param_bounds(i,2)-param_bounds(i,1)) * lhs_unit(:, i);
end

%% Load in the observed formation & surface age distributions
% Formation age distribution data from Korenaga (2018a)
%data_formationage = load('korenaga18a_Tunmix_orig.dat');
data_formationage = readmatrix(...
    'cc/observation/cc_age.xlsx','Sheet','korenaga18a_Tunmix_orig');
% Surface age distribution data from Roberts & Spencer (2015)
data_zircon_surf = readmatrix(...
    'cc/observation/cc_age.xlsx','Sheet','korenaga18a_T_U_Pb');
%[F_Jun_same] = load_FandS_fun(t,nt,data_formationage);
% use load_FandS_fun function to set data in the same dimension as the time series
[F_Jun_same,S_Jun_same] = load_FandS_fun(t,nt,data_formationage,data_zircon_surf);


cluster_end=floor(N_all/N_one_time);

futures(N_one_time) = parallel.FevalFuture;
results = zeros(N_all,8);

for cluster=1:cluster_end
    Ns=1+(cluster-1)*N_one_time;
    Ne=cluster*N_one_time;
    % Nl=Ne-Ns+1;   
    k=0;
    for i = Ns:Ne
            k=k+1;
            disp(i);
           % disp(max(t));
            param_struct = struct('kappa_r', params(i,1), ...
                          'kappa_g', params(i,2), ...
                          'Rs',      params(i,3), ...
                          'Rp',      params(i,4), ...
                          'ts',      params(i,5), ...
                          'frw',     params(i,6), ...
                          'Mcp',     Mcp, ...
                          't',       t, ...
                          'tmax',    tmax, ...
                          'F_obs',   F_Jun_same, ...
                          'S_obs',   S_Jun_same);
    futures(k) = parfeval(@run_cc_growth_model, 1, param_struct);        % end
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
            results(i, :) = nan(1,8);  % 
        end
    end


end


elapsedTime = toc;  % s
disp(['used time = ', num2str(elapsedTime), ' s']);

% results(:,1)=results(:,1)/1e9/yr_s; % time_cc from s to Ga
col_names = {'ts','kappa_g', 'kappa_r', 'Rs', 'Rp', 'frw', 'RMSE_F', 'RMSE_S'};

save('cc/out_cc_growth_lhs_Krw_lowf.mat', 'results', 'col_names','-v7.3');

%% postprocesses
%% select successful solutions for plotting results
crit_misfit_F = 0.1;
crit_misfit_S = 0.2;

rmseF_idx = find(strcmp(col_names, 'RMSE_F'));
rmseS_idx = find(strcmp(col_names, 'RMSE_S'));
success_mask = results(:, rmseF_idx) < crit_misfit_F & results(:, rmseS_idx) < crit_misfit_S;
success_F = results(success_mask, :);

% whether sampling number is enough
if size(success_F, 1) < 1000
    disp('Success sample not enough');
end

save('cc/out_cc_grwoth_success_lhs_MC_Krw_lowf.mat', 'success_F', '-v7.3');
