

%% prepare data for misfit plane:
% tune [kgcc or Rs or Rp] fixed other parameters to get countours with FrXe
% obtained data for misfit plane plot
% clear;clc;

%% load successful seeds from stage1-3 

success_F_struct=load('cc/out_cc_grwoth_success_lhs_MC_Krw_lowf.mat');
success_F=success_F_struct.success_F;
[size_F,nparam_cc]=size(success_F);

success_T_struct=load('thermal/out_thermal_success_lhs_test_fixK_U_MC_lowerT_Rac1100_lowf.mat');
success_T=success_T_struct.success_T;
[size_T,nparam_T]=size(success_T);

success_Xe_struct=load('Xe/out_Xe_success_lhs_new_3e7_test_fixK_U_MC_negativFD_Rac1100_new_lowf.mat');
success_Xe=success_Xe_struct.success_Xe;
[size_Xe,nparam_Xe]=size(success_Xe);

%%
% j=randi(size_Xe);
j=2000;

% %%
% % if not in successful
% % BSE_set=-3;
% % FrXe=0:0.05:0.5;
% % idbse=find(seeds_stage4(:,6)==BSE_set);
% % j=idbse(1);
% 
% BSE_set=-1;
% 
% seeds_c_struct=load('water/out_for_seeds_stage4.mat');
% seeds_c=seeds_c_struct.seeds_stage4_sorted;
% seeds_stage4 = seeds_c;
% 
% idbse1=find(seeds_c(:,6)==BSE_set);
% j=idbse1(1);%for this bse, fixed T curve

%% set params
par_solver_new;
tmax = t_pd/yr_s/1e9;% the age of solar system, in unit Ga

Mcp = mcc_val;% mass of continental crust at present-day, in unit kg
para_Scale=2;% 2: Schubert
Msi=1.5;% 
SwXe = logical(1);
Xei=3.2e7; % atoms/g

%% set input
default_input_for_degassing;
t=Simulation_input.time_series;
t=t/yr_s/1e9;% to Ga
nt = length(t);
ntd=nt;
td=t;

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


%% set FrXe & cc variables
FrXe=0:0.05:0.5;

kappa_gcc=success_Xe(j,2);
kappa_rcc=success_Xe(j,3);
Rs=success_Xe(j,4);
Rp=success_Xe(j,5);
Ubse=success_Xe(j,6);

% kappa_gcc = -1:3.2:31;% cannot have 0, kgcc=0, error=nan (Mc eq has 1/0)
% Yq=kappa_gcc;

Ubse = 10e-9:1e-9:20e-9;
Yq=Ubse;

% Rs= 0:2e22:20e22;
% Yq=Rs;

% kappa_rcc = -3:0.6:3;% cannot have 0, kgcc=0, error=nan (Mc eq has 1/0)
% Yq=kappa_rcc;

% Rp= 0:0.2e22:2e22;
% Yq=Rp;

% 
% kappa_gcc = 30;
%  Rp=1.2e22;
%  Rs=6e22;


%%
ncores=6;
if isempty(gcp('nocreate'))
    parpool('local', ncores);  % 
end

%% Submit parfeval jobs
%futures(length(FrXe)*50) = parallel.FevalFuture;
futures(length(FrXe)*length(Yq)) = parallel.FevalFuture;
k=0;

tic;
for i = 1:length(FrXe)
    for n = 1:length(Yq)
        k=k+1;
        disp(k);
       param_struct = struct( 't',       t*1e9*yr_s, ... % change from Ga to s
                              'tmax',    tmax*1e9*yr_s, ... % change from Ga to s ====== above: time vector =======
                              'time_cc_input', success_Xe(j,1)*1e9*yr_s,... %seeds_stage4(j,1)*1e9*yr_s,... % 0.2*1e9*yr_s, ... %   0.2*1e9*yr_s, ... % change from Ga to s
                              'kappa_gcc', kappa_gcc,... % success_Xe(j,2),...% seeds_stage4(j,2), ... %  6,... % 
                              'kappa_rcc', kappa_rcc,... %success_Xe(j,3),...%seeds_stage4(j,3), ... %  1,... %
                              'Rs',  Rs,... %success_Xe(j,4),... % seeds_stage4(j,4), ... %  10e22, ... % 
                              'Rp', success_Xe(j,5),...%seeds_stage4(j,5), ... %  1.5e22, ... % % ====== above: CC model =======
                              'BSE',   -4,...  % seeds_stage4(j,6), ...%  4,...% 
                              'Ubse', Ubse(n),... %success_Xe(j,6),...%seeds_stage4(j,6), ...% 
                              'Qc_pd_input',       success_Xe(j,7),...%seeds_stage4(j,7), ... %  15e12, ... %
                              'dQc_input',      0,...%seeds_stage4(j,8), ... % 5e12,... % 
                              'Ti',        success_Xe(j,9),...% seeds_stage4(j,9),...    % 2500,...%                          
                              'eta_ref',        success_Xe(j,10),...%seeds_stage4(j,10), ... %  7.5e14,...% % ====== above: Thermal model ====== 
                              'Fd_mor',    success_Xe(j,16),...  % 0.5,...%seeds_stage4(j,16),... %1, ... % 
                              'Fd_p',     success_Xe(j,17),...% 0.5,...%%seeds_stage4(j,17),... % 1, ... %  % ====== above: MOR & Plume degassing ====== 
                              'P_mor_pd', P_mor_pd,...
                              'Fr_w',     success_Xe(j,13),...% seeds_stage4(j,13),...  % 0.25, ...       
                              'Xm_init',  success_Xe(j,11),...% seeds_stage4(j,11),... %  0, ... % 
                              'Ms_init',    success_Xe(j,12),...% seeds_stage4(j,12),...
                              'water_obs1', water_obs1,... % ====== above: watre model ====== 
                              'SwXe',  SwXe,...
                              'Fr_Xe',      FrXe(i),...                             
                              'Xe_init',    success_Xe(j,14),...                             
                              'Xe_obs', Xe_obs ... % ====== above: Xe model ======
                               );

        futures(k) = parfeval(@run_degassing_model, 1, param_struct);
        %out_test=run_degassing_model(param_struct);
    end
end


%% Collect outputs
results = zeros(length(FrXe)*length(Yq),25);
for i = 1:length(FrXe)*length(Yq)%50
    if mod(i,10)==0
        disp(['lccFr=' num2str(i) ' of ' num2str(length(FrXe)*length(Yq))]);% keep track of the calculation
    end
    try
        results(i, :) = fetchOutputs(futures(i));
    catch ME
        warning('Error in task %d: %s', i, ME.message);
        results(i, :) = nan(1,25);  % 
    end
end

% results = zeros(length(FrXe)*length(Yq),25);
% for i = 1:length(FrXe)*length(Yq)
%     try
%         results(i, :) = fetchOutputs(futures(i));
%     catch ME
%         warning('Error in task %d: %s', i, ME.message);
%         results(i, :) = nan(1,25);  % 
%     end
% end

%% Save results
results(:,1)=results(:,1)/1e9/yr_s; % time_cc from s to Ga
col_names = {'ts_cc','kappa_gcc','kappa_rcc','Rs','Rp',...
    'Ubse','Qc_pd', 'dQc', 'Ti', 'eta_ref',...
    'Xm_init','Ms_init','Fr_w',...
    'Xe_init','Fr_Xe',...
    'Fd_mor','Fd_p',...
    'RMSE_water',...
    'Er_xec','Er_xe128v130','Er_xe128v132','Er_xe130v132','Er_xe131v132','Er_xe134v132','Er_xe136v132'};
%    'RMSE_xec','RMSE_xe128v130','RMSE_xe128v132','RMSE_xe130v132','RMSE_xe131v132','RMSE_xe134v132','RMSE_xe136v132'};

% Ubse_Fr_Xe_contours=results;
% krcc_Fr_Xe_contours=results;
% kgcc_Fr_Xe_contours=results;
% Rs_Fr_Xe_contours=results;

% kappa_gcc = -1:3.2:31;% cannot have 0, kgcc=0, error=nan (Mc eq has 1/0)
% 
% 
% Ubse = 10e-9:1e-9:20e-9;
% 
% 
% Rs= 0:2e22:20e22;
% 
% kappa_rcc = -3:0.6:3;% cannot have 0, krcc=0, error=nan (Mc eq has 1/0)


% save('Xe/fin_new_2000_RpFr_Xe_3e7.mat', 'kappa_gcc','kappa_rcc','Rs','Ubse','Ubse_Fr_Xe_contours','kgcc_Fr_Xe_contours', 'col_names','FrXe','-v7.3');
% save('Xe/fin_new_2000_RpFr_Xe_3e7.mat', 'kappa_gcc','kappa_rcc','Rs','Ubse','Ubse_Fr_Xe_contours','krcc_Fr_Xe_contours','kgcc_Fr_Xe_contours','Rs_Fr_Xe_contours', 'col_names','FrXe','-v7.3');

% delete(gcp('nocreate'));

%% postprocesses
%% select successful solutions for plotting results
crit_misfit_w = 1;
crit_misfit_xec = 1;
crit_misfit_xe128v130 = 1; 
crit_misfit_xe128v132 = 1; 
crit_misfit_xe130v132 = 1; 
crit_misfit_xe131v132 = 1; 
crit_misfit_xe134v132 = 1; 
crit_misfit_xe136v132 = 1; 

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
                    abs( results(:, Erxe136v132_idx) ) <= crit_misfit_xe136v132 ); %& ...
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

% success_Xe = results(success_mask, :);
% 
% % whether sampling number is enough
% if size(success_Xe, 1) < 1
%     disp('No Success Sample');
% end

%% === overlap ===
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

subset1 = results;

Erxec_e = subset1(:, Erxec_idx);
Erxe128v130_e = subset1(:, Erxe128v130_idx);
Erxe130v132_e = subset1(:, Erxe130v132_idx);
Erxe131v132_e = subset1(:, Erxe131v132_idx);
Erxe134v132_e = subset1(:, Erxe134v132_idx);
Erxe136v132_e = subset1(:, Erxe136v132_idx);

 Xq=FrXe;

Erxec_e_grid = reshape(Erxec_e,[11,11]);
Erxe128v130_e_grid = reshape(Erxe128v130_e,[11,11]);
Erxe130v132_e_grid = reshape(Erxe130v132_e,[11,11]);
Erxe131v132_e_grid = reshape(Erxe131v132_e,[11,11]);
Erxe134v132_e_grid = reshape(Erxe134v132_e,[11,11]);
Erxe136v132_e_grid = reshape(Erxe136v132_e,[11,11]);

% [Xq, Yq] = meshgrid(FrXe, Yq);
%  X=FrXe;
%  Y=Yq;
% Erxec_e_grid = griddata(X, Y, Erxec_e, Xq, Yq, 'natural');
% Erxe128v130_e_grid = griddata(X, Y, Erxe128v130_e, Xq, Yq, 'natural');
% Erxe130v132_e_grid = griddata(X, Y, Erxe130v132_e, Xq, Yq, 'natural');
% Erxe131v132_e_grid = griddata(X, Y, Erxe131v132_e, Xq, Yq, 'natural');
% Erxe134v132_e_grid = griddata(X, Y, Erxe134v132_e, Xq, Yq, 'natural');
% Erxe136v132_e_grid = griddata(X, Y, Erxe136v132_e, Xq, Yq, 'natural');
% 

%% plot countour
% [Xpd Ypd] = contourf(...
%     FrXe,Yq,...
%     Erxec_e_grid);colorbar;
% 
% [Xpd Ypd] = contourf(...
%     FrXe,Yq,...
%    Erxe128v130_e_grid);colorbar;

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

figure;
set(gcf,'color','w');
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


%%

elapsedTime = toc;  % s
disp(['used time = ', num2str(elapsedTime), ' s']);





