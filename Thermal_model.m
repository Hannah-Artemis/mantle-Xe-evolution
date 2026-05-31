function [Tp,Tm,Qs,Qc,H,Ra,eta] = Thermal_model(t,Ti,eta_ref,Qc_pd_input,dQc_input,...
    BSE,Ubse,kappa_gcc,time_cc_input,para_Scale,Hsf)

%% % initialize
par_solver_new;

t = t; % in the unit [s]
dt = t(2)-t(1);% calcualte timestep dt
nt=length(t);
Tm = nan(nt,1); % average mantle temperature
Tp = nan(nt,1); % average mantle temperature
Qs = nan(nt,1); % surface heat flux
Ra = nan(nt,1); 
H = nan(nt,1); 
H_El = nan(nt,4);

R0U238_U = 0.9927;
R0U235_U = 1-R0U238_U;
R0Th_U = 4;% 4[GK-CI]; 3.08; 3.0-3.8[2014-Javoy-EN]
R0K_U = 1.27e4;% 1.27e4[GK];1.21e4
R0K40_K = 1.28e-4;% 1.28e-4[GK];1.19e-4


U238bse = Ubse*R0U238_U;
U235bse = Ubse*R0U235_U;
Thbse =   Ubse*R0Th_U;
K40bse =  Ubse*R0K_U*R0K40_K;


h1 = [U238bse;  9.46e-5;    4.4679e9*yr_s     ]; %U238 C30.8 MS13 E1-10.7 E2-6.31
h2 = [U235bse;  5.69e-4;    0.704e9*yr_s    ];  %U235 0.22 0.09 0.0757 0.0445
h3 = [Thbse;   2.64e-5;    14e9*yr_s       ];    %Th 124 48.6 23.1 10.8
h4 = [K40bse;  2.92e-5;    1.25e9*yr_s     ]; %K 36.9 18.7 55.6 36
h_El = [h1,h2,h3,h4];




cc = (t<time_cc_input).*0+(t>=time_cc_input)...
    .*(1-exp(-kappa_gcc*(t-time_cc_input)/yr_s/1e9))...
    ./(1-exp(-kappa_gcc*(t_pd-time_cc_input)/yr_s/1e9 ) ) ;


Hsf_t=1;
H_El(:,1) = (U238bse-U238CC*mcc_val*(cc)/rho/V/Hsf_t).*h_El(2,1).*exp((log(2).*(t_pd-t))./h_El(3,1));
H_El(:,2) = (U235bse-U235CC*mcc_val*(cc)/rho/V/Hsf_t).*h_El(2,2).*exp((log(2).*(t_pd-t))./h_El(3,2));
H_El(:,3) = (Thbse-ThCC*mcc_val*(cc)/rho/V/Hsf_t).*h_El(2,3).*exp((log(2).*(t_pd-t))./h_El(3,3));
H_El(:,4) = (K40bse-KCC*mcc_val*(cc)/rho/V/Hsf_t).*h_El(2,4).*exp((log(2).*(t_pd-t))./h_El(3,4));

% use Qc_pd to calculate hsf
H = sum(H_El,2)*rho*V;
Hsf = 1+Qc_pd_input/H(end);
Qc = sum(H_El,2)*(Hsf-1)*rho*V;

%% % initial value

Tm(1)=Ti;
Ra(1) = (alpha_ra.*rhoa.*g.*(Tm(1)-Ts).*d^3)./(k.*(eta_ref.*exp(E./(Rg.*Tm(1)))));
switch para_Scale
    case 1 % Crowley, 2011: using weakly nonlinear theory from Geodynamics      
        Qs(1) = 2*(S+(1-cc(1))*Scc)*Kc.*Tm(1).*(U./(pi*L*k)).^(1/2); %W
    case 2 % Schubert, using empirical critical Ra (not for onset but for stable convection)
        Qs(1) = 2*(S+(1-cc(1))*Scc)*Kc.*(Tm(1)-Ts)/d.*(Ra(1)/Ra_cr).^(1/3); %   
    case 3 % Labrosse, 2007: using fixed max plate age & triangle age distribution                  
        Qs(1) = (2).*(S)*Kc.*Tm(1).*(1./(pi*age_max*k)).^(1/2); %W
    case 4 % Labrosse, 2007: using fixed max plate age 
        %Qs = (2+cc*2/3).*(S+(1-cc)*Scc)*Kc.*T.*(1./(pi*age_max*k)).^(1/2); %W
        Qs(1) = (8/3).*(S)*Kc.*Tm(1).*(1./(pi*age_max*k)).^(1/2); %W
end



% convert average temperature Tm [K] to potential temperature Tp [C]
Tp = Tm.*exp(-alpha_ra*g*1445e3/Cp)-273;
% Tp(1) = Tm(1).*exp(-alpha_ra*g*1445e3/Cp)-273;
eta=eta_ref.*exp(E./(Rg.*Tm));

%% % evolution
for i = 2:nt
    Tp(i) = Tp(i-1)+((-Qs(i-1)+H(i-1)+Qc(i-1))/(Cp*rho*V))*(t(i)-t(i-1));
    Tm(i) = (Tp(i)+273)./exp(-alpha_ra*g*1445e3/Cp);
    Ra(i) = (alpha_ra.*rhoa.*g.*(Tm(i)-Ts).*d^3)./(k.*(eta_ref*exp(E/Rg/Tm(i))));
    switch para_Scale
        case 1 % Crowley, 2011: using weakly nonlinear theory from Geodynamics      
            Qs(i) = 2*(S+(1-cc(i))*Scc)*Kc.*Tm(i).*(U./(pi*L*k)).^(1/2); %W
        case 2 % Schubert, using empirical critical Ra (not for onset but for stable convection)
            Qs(i) = 2*(S+(1-cc(i))*Scc)*Kc.*(Tm(i)-Ts)/d.*(Ra(i)/Ra_cr).^(1/3); %   
        case 3 % Labrosse, 2007: using fixed max plate age & triangle age distribution                  
            Qs(i) = (2).*(S)*Kc.*Tm(i).*(1./(pi*age_max*k)).^(1/2); %W
        case 4 % Labrosse, 2007: using fixed max plate age 
            %Qs = (2+cc*2/3).*(S+(1-cc)*Scc)*Kc.*T.*(1./(pi*age_max*k)).^(1/2); %W
            Qs(i) = (8/3).*(S)*Kc.*Tm(i).*(1./(pi*age_max*k)).^(1/2); %W
    end
end

end