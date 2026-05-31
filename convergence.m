%% % convergence test


%%
default_input_for_degassing;
success_Xe_struct=load('Xe/out_Xe_success_lhs_new_3e7_test_fixK_U.mat');
success_Xe=success_Xe_struct.success_Xe;

% [size_Xe,nparameter]= size(success_Xe);

size_Xe = 1;
% success_Xe = success_Xe(4:end,:);

l=1;

param_success = struct ('time_cc_input', success_Xe(l,1)*1e9*yr_s,...% here time_cc is in s
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
                            'Fr_w', success_Xe(l,13),...
                            'Xe_init', success_Xe(l,14),...
                            'Fr_Xe', success_Xe(l,15),...
                            'Fd_mor', success_Xe(l,16),...
                            'Fd_p', success_Xe(l,17),...
                            't',       t*1e9*yr_s, ... % change from Ga to s
                            'tmax',    t_pd, ... % change from Ga to s ====== above: time vector =======
                            ... 'water_obs1', water_obs1,... % ====== above: watre model ====== 
                            'SwXe',  logical(1) ...
                            ...'Xe_obs', Xe_obs ... % ====== above: Xe model ======
                             );

[CC_growth_input,Thermal_input,H2O_input,Xe_input,Simulation_input,FIG,Fdmor_model,Fdp_model] = input_degassing(param_success);



 

%%
ntd=460001;

Simulation_input.nsteps  = Simulation_input.end*[ ntd ]; % number of timesteps
    t_fin = Simulation_input.end*t_pd;
    t = linspace(0,t_fin,Simulation_input.nsteps);
    Simulation_input.time_series = t';

NGcc = nan(ntd,size_Xe);
Rcc = nan(ntd,size_Xe);
Gcc = nan(ntd,size_Xe);

Tp = nan(ntd,size_Xe);
Qs = nan(ntd,size_Xe);

% Xm_water = nan(ntd,size_Xe);

Xe = nan(ntd,size_Xe);
Xe128 = nan(ntd,size_Xe);
Xe131 = nan(ntd,size_Xe);
Xe132 = nan(ntd,size_Xe);
Xe134 = nan(ntd,size_Xe);
Xe136 = nan(ntd,size_Xe);

Pu = nan(ntd,size_Xe);
Ur = nan(ntd,size_Xe);

 [NGcc_model,Rcc_model,Gcc_model] = CC_growth_fun1(t,CC_growth_input.time_cc/1e9/yr_s,t_pd/1e9/yr_s,...
        mcc_val,CC_growth_input.kcc,CC_growth_input.Rccp,CC_growth_input.Rccs,CC_growth_input.krcc);% here time_cc is in Ga
        % save the results
        NGcc(:,l) = NGcc_model;
        Rcc(:,l) = Rcc_model;
        Gcc(:,l) = Gcc_model;
      
% [Tp_model,Tm_model,Qs_model,Qc_model,H_model,Ra_model,eta_out] = Thermal_model(td*1e9*yr_s,param_success.Ti,param_success.eta_ref,param_success.Qc_pd_input,param_success.dQc_input,...
%     param_success.BSE,param_success.Ubse,param_success.kappa_gcc,param_success.time_cc_input,para_Scale,Thermal_input.Hsf);

 output_model = Degassing_model(Simulation_input,...
                     Thermal_input,CC_growth_input,H2O_input,Xe_input,Fdmor_model,Fdp_model);
        
% save the results
Tp(:,l) = output_model.T*exp(-alpha_ra*g*1445e3/Cp)-273; %C
Qs(:,l) = output_model.Qs;
%        Xm_water(:,l) =  output_model.Xm;
Xe(:,l) =  output_model.Xe;
Xe128(:,l) = output_model.Xe_Atm;
Xe131(:,l) = output_model.Xe131;
Xe132(:,l) = output_model.Xe132;
Xe134(:,l) = output_model.Xe134;
Xe136(:,l) = output_model.Xe136;
Pu(:,l) = output_model.Put;
Ur(:,l) = output_model.Urt;

%%
ntd=4601;

Simulation_input.nsteps  = Simulation_input.end*[ ntd ]; % number of timesteps
    t_fin = Simulation_input.end*t_pd;
    t = linspace(0,t_fin,Simulation_input.nsteps);
    t = [0:0.0001:2,2.001:0.001:2.6,2.605:0.005:4.6]/4.6*t_pd;
    ntd=length(t);
    Simulation_input.time_series = t';


NGcc1 = nan(ntd,size_Xe);
Rcc1 = nan(ntd,size_Xe);
Gcc1 = nan(ntd,size_Xe);

Tp1 = nan(ntd,size_Xe);
Qs1 = nan(ntd,size_Xe);

% Xm_water = nan(ntd,size_Xe);

Xe1 = nan(ntd,size_Xe);
Xe1281 = nan(ntd,size_Xe);
Xe1311 = nan(ntd,size_Xe);
Xe1321 = nan(ntd,size_Xe);
Xe1341 = nan(ntd,size_Xe);
Xe1361 = nan(ntd,size_Xe);

Pu1 = nan(ntd,size_Xe);
Ur1 = nan(ntd,size_Xe);


[NGcc_model,Rcc_model,Gcc_model] = CC_growth_fun1(t,CC_growth_input.time_cc/1e9/yr_s,t_pd/1e9/yr_s,...
        mcc_val,CC_growth_input.kcc,CC_growth_input.Rccp,CC_growth_input.Rccs,CC_growth_input.krcc);% here time_cc is in Ga
        % save the results
        NGcc1(:,l) = NGcc_model;
        Rcc1(:,l) = Rcc_model;
        Gcc1(:,l) = Gcc_model;

output_model1 = Degassing_model(Simulation_input,...
             Thermal_input,CC_growth_input,H2O_input,Xe_input,Fdmor_model,Fdp_model);
        
% save the results
Tp1(:,l) = output_model1.T*exp(-alpha_ra*g*1445e3/Cp)-273; %C
Qs1(:,l) = output_model1.Qs;
%        Xm_water(:,l) =  output_model1.Xm;
Xe1(:,l) =  output_model1.Xe;
Xe1281(:,l) = output_model1.Xe_Atm;
Xe1311(:,l) = output_model1.Xe131;
Xe1321(:,l) = output_model1.Xe132;
Xe1341(:,l) = output_model1.Xe134;
Xe1361(:,l) = output_model1.Xe136;
Pu1(:,l) = output_model1.Put;
Ur1(:,l) = output_model1.Urt;
   
%% % plot
% figure(111);
% set(gcf,'color','white');
% subplot(2,1,1); hold off;
% plot(td,(Xe-Xe1)./(Xe+Xe1)*2,'k-');
% xlabel('Time (Gyr)'); ylabel('relative error of Xe'); grid on;
% subplot(2,1,2); hold off;
% plot(td,(Xe1361-Xe1361)./(Xe136+Xe1361)*2,'k-');
% xlabel('Time (Gyr)'); ylabel('relative error of 136Xe'); grid on;
% subplot(3,1,3); hold off;
% plot(td,(Pu-Pu1)./(Pu+Pu1)*2,'k-');
% xlabel('Time (Gyr)'); ylabel('relative error of Xe'); grid on;

% figure(115);
% set(gcf,'color','white');
% plot(t,Ur_50);


disp('cc');
disp((NGcc(end)-NGcc1(end))./(NGcc(end)+NGcc1(end))*2);
disp((Gcc(end)-Gcc1(end))./(Gcc(end)+Gcc1(end))*2);
disp('T');
disp((Tp(end)-Tp1(end))./(Tp(end)+Tp1(end))*2);
disp((Qs(end)-Qs1(end))./(Qs(end)+Qs1(end))*2);
disp('Xe');
disp((Xe(end)-Xe1(end))./(Xe(end)+Xe1(end))*2);
disp((Xe136(end)-Xe1361(end))./(Xe136(end)+Xe1361(end))*2);




