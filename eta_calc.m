function eta = eta_calc(Xm,T,etalawin,etain)
   
% if(nargin < 5)
%    etain = 8.73e14;
% end

par_solver_new;
global f_WI

 Rf =    - [0 0.004605 0.00921 0.01382 0.01842 0.02303];
%Rf = -[0 0.00115 0.0230];

switch etalawin
    % papers
    case 0 % CONSTANT
        eta = etain.*exp(0.*T);
        
    case 1  %'McG'
        eta=8.9e10.*exp((McG_E)./(Rg.*T)); % From McGovern and Schubert (1980)
        
    case 2  %'AnB'
        Xm_wtfr = Xm*1e-6; % ppm --> weight fraction
        eta=8.9e10*exp((6.4e4 + (-6.1e6.*Xm_wtfr))./T); % McG1989 Anita Bay parametrisation
        
    case 3  %'Ahm'
        Xm_wtfr = Xm*1e-6; % ppm --> weight fraction
        eta=8.9e10*exp((6.4e4 + (-8.1e5.*Xm_wtfr))./T);  % McG1989 Aheim parametrisation
        
    case 4  %'Cro'
        Coh = Xm * 100/6; % ppm --> H/10^6Si
        F_h2o = exp(c0 + (c1*log(Coh)) + (c2*log(Coh).*log(Coh)) + (c3*log(Coh).*log(Coh).*log(Coh))); %MPa to Pa
        pe = A*(F_h2o.^(-r_eta));% (Pas/MPa-r)MPa^-r
        eta = pe.* exp(E./(Rg.*T)); %Pas
        
    case 5  %'San'
        Coh = Xm * 100/6; % ppm --> H/10^6Si
        F_h2o = exp(c0 + (c1*log(Coh)) + (c2*log(Coh).*log(Coh)) + (c3*log(Coh).*log(Coh).*log(Coh))); %MPa to Pa
        pe = San_eta0/San_Acre*(F_h2o.^(-San_r_eta));% 
        eta = pe.* exp(San_E./(Rg.*T));
        
    case 6 % Cro uncorrected
        Coh = Xm * 6.3e-4; % ppm --> H/10^6Si
        F_h2o = 1e6*exp(c0 + (c1*log(Coh)) + (c2*log(Coh).*log(Coh)) + (c3*log(Coh).*log(Coh).*log(Coh))); %MPa to Pa
        pe = A_un*F_h2o.^-(r_eta/n_eta);
        eta = pe.* exp(E./(Rg.*T));
        
    case 7 % Cro uncorrected dry case
        Coh = f_WI * 6.3e-4; % 300 ppm --> H/10^6Si
        F_h2o = 1e6*exp(c0 + (c1*log(Coh)) + (c2*log(Coh).*log(Coh)) + (c3*log(Coh).*log(Coh).*log(Coh))); %MPa to Pa
        pe = A_un*F_h2o.^-(r_eta/n_eta);
        eta = pe.* exp(E./(Rg.*T));
        
    case 8 % Crowley code T corrected WI
        Coh = f_WI * 6.3e-4; % 300 ppm --> H/10^6Si
        F_h2o = 1e6*exp(c0 + (c1*log(Coh)) + (c2*log(Coh).*log(Coh)) + (c3*log(Coh).*log(Coh).*log(Coh))); %MPa to Pa
        pe = A_un*F_h2o.^-(r_eta/n_eta);
        eta = pe.* exp(E./(Rg.*T));    
        
    case 9 % Cro correct dry case
        Coh = f_WI * 100/6; % ppm --> H/10^6Si
        F_h2o = exp(c0 + (c1*log(Coh)) + (c2*log(Coh).*log(Coh)) + (c3*log(Coh).*log(Coh).*log(Coh))); %MPa to Pa
        pe = A*(F_h2o.^(-r_eta));% (Pas/MPa-r)MPa^-r
        eta = pe.* exp(E./(Rg.*T)); %Pas
        
        
        
         % own laws
         % *** 500 ppm calibration ***
    case {10 1000}
        eta = etain.*exp(Rf(1).*Xm).*exp(E./(Rg.*T)); %8.73e14
        
    case 12
        eta = etain.*exp(Rf(2).*Xm).*exp(E./(Rg.*T)); %8.73e15
        
    case 14
        eta = etain.*exp(Rf(3).*Xm).*exp(E./(Rg.*T)); %8.73e16
        
    case 15
        eta = etain.*exp(-0.0106.*Xm).*exp(E./(Rg.*T)); %7.53e16
        
    case 16
        eta = etain.*exp(Rf(4).*Xm).*exp(E./(Rg.*T)); %8.73e17
        
    case 18
        eta = etain.*exp(Rf(5).*Xm).*exp(E./(Rg.*T)); % 8.73e18
        
    case 20
        eta = etain.*exp(Rf(6).*Xm).*exp(E./(Rg.*T)); %8.73e19
        
    case 21
        eta = etain.*exp(-0.02303.*Xm).*exp(E./(Rg.*T)); %8.73e16
        
        % *********************  200 ppm calibration ********************* 
    case 200 % O(0) 
        eta = etain .* exp(0.*Xm).*exp(E./(Rg.*T));

    case 202 % O(~2)
        eta = 3.49e15 .* [exp(-0.0139.*Xm) + (0.1054.*exp(-0.0019.*Xm))] .* exp(E./(Rg.*T));
    
    case 204 % O(~4)
        eta = 3.77e17 .* [exp(-0.0346.*Xm) + (0.0003.*exp(-0.0023.*Xm))] .* exp(E./(Rg.*T));
        
       % ********************* 300 ppm calibration ********************* 
    case 302 % O(~2)
        eta = 3.54e15 .* [exp(-0.0087.*Xm) + (0.0785.*exp(-0.0016.*Xm))] .* exp(E./(Rg.*T));
        
    case 304 % O(~4)
        eta = 6.25e13 .* [exp(-0.0015.*Xm) + (6022.5.*exp(-0.0230.*Xm))] .* exp(E./(Rg.*T));
            
    case 306 % O(~6)
        eta = 3.77e19 .* [exp(-0.0384.*Xm) + (8.3049e-7.*exp(-0.0014.*Xm))] .* exp(E./(Rg.*T));
        
       % ********************* 400 ppm calibration ********************* 
    case 402 % O(~2)
        eta = 3.67e15 .* [exp(-0.0061.*Xm) + (0.0321.*exp(-0.0010.*Xm))] .* exp(E./(Rg.*T));
        
    case 404 % O(~4)
        eta = 2.31e13 .* [exp(-8.1108e-4.*Xm) + (16321.*exp(-0.0173.*Xm))] .* exp(E./(Rg.*T));
            
    case 407 % O(~7)
        eta = 3.77e19 .* [exp(-0.0288.*Xm) + (3.0101e-7.*exp(-7.0562e-4.*Xm))] .* exp(E./(Rg.*T));
        
        % *********************  150 ppm calibration ********************* 
    case 1500 % O(0) 
        eta = etain .* exp(0.*Xm).*exp(E./(Rg.*T));

    case 1502 % O(~2)
        eta = 3.49e15 .* [exp(-0.0191.*Xm) + (0.1101.*exp(-0.0020.*Xm))] .* exp(E./(Rg.*T));
    
    case 1504 % O(~4)
        eta = 3.77e16 .* [exp(-0.0313.*Xm) + (0.0065.*exp(-0.0026.*Xm))] .* exp(E./(Rg.*T));
        
        % *********************  1e23 0ppm Calibration ********************* 
        % Where laws 'matching' experimental data
        
    case 100 % no water influence
        eta = etain .* exp(E./(Rg.*T));
    case 99 % no water influence
        eta = etain .* exp(E./(Rg.*T));
    case 101 % simple straight line O(1)
        eta = 7.53e+15 .* exp(-0.0018.*Xm) .* exp(E./(Rg.*T));
    case 103 % simple straight line O(3)
        eta = 7.53e+15 .* exp(-0.0067.*Xm) .* exp(E./(Rg.*T));
    case 121 % power law O(1)
        eta = etain .* (Xm./1).^-0.3 .* exp(E./(Rg.*T));
    case 123 % power law O(3)
        eta = etain .* (Xm./1).^-1.0 .* exp(E./(Rg.*T));
        
        % *********************  5e21 456ppm Calibration ********************* 
        % Where laws to match 5e21 average mantle viscosity (geoid)
        
    case 1001 % simple straight line O(1)
        eta = 7.63e14  .* exp(-0.0018.*Xm) .* exp(E./(Rg.*T));
    case 1003 % simple straight line O(3)
        eta = 7.53e15 .* exp(-0.0066.*Xm) .* exp(E./(Rg.*T));
    case 1201 % power law O(1)
        eta = 3.77e14 .* (Xm./500).^-0.3 .* exp(E./(Rg.*T));
    case 1203 % power law O(3)
%         Xm = max(Xm(:,1),1);
        eta = 3.77e14 .* (Xm./500).^-1.0 .* exp(E./(Rg.*T));
        
        % *********************  1e23 500ppm Calibration ********************* 
    case 2350
        eta = etain .* exp(E./(Rg.*T));
    case 2351
        eta = 7.53e15 .* (Xm./500).^-0.3 .* exp(E./(Rg.*T));
    case 2353
        eta = 7.53e15 .* (Xm./500).^-1 .* exp(E./(Rg.*T));
    case 2350400
        eta = etain .* exp(400e3./(Rg.*T));
    case 2350500
        eta = etain .* exp(500e3./(Rg.*T));
        
        
        
%     case 1001 % Nakagawa Law weak
%         eta = 3.765e14 .* (Xm./1).^-0.3 .* exp(E./(Rg.*T)); 
%         
%     case 1003 % Nakagawa Law strong
%         eta = 3.765e14 .* (Xm./1).^-1 .* exp(E./(Rg.*T));
%         
%     case 2001 % Simple Law O(1) 500ppm calibration
% %         eta = 7.6436e+14 .* exp(-0.0017.*Xm) .* exp(E./(Rg.*T));
%         eta = 7.63e+14 .* exp(-0.0018.*Xm) .* exp(E./(Rg.*T));
%     case 2003 % Simple Law O(3) 500ppm calibration
% %         eta = 7.5307e+15 .* exp(-0.006.*Xm)  .* exp(E./(Rg*T));
%         eta = 7.53e+15 .* exp(-0.0067.*Xm)  .* exp(E./(Rg*T));
        
%     otherwise 
%         error(['Unknown etalaw: ' etalawin])
end

%     case 10 % Our own code O(11) change over 0-1000ppm
%         eta = 7.5e18.*exp(-0.02533.*Xm).*exp(E./(Rg.*T)); %7.5
%         
%     case 11 % Our own code O(2) change over 0-1000ppm
%         eta = etain.*1e1.*exp(-0.004605.*Xm).*exp(E./(Rg.*T)); %8.73e15
%         
%     case 12 % Our own code O(1) change over 0-1000ppm
%         eta = 3.75e15.*exp(-0.002303.*Xm).*exp(E./(Rg.*T));
%         
%     case 13 % Our own code O(8) change over 0-1000ppm
%         eta = etain.*1e4.*exp(-0.01842.*Xm).*exp(E./(Rg.*T)); % 8.73e18
%         
%     case 14 % Our own code O(6) change over 0-1000ppm
%         eta = etain.*1e3.*exp(-0.01382.*Xm).*exp(E./(Rg.*T)); %8.73e17
%         
%     case 15 % Our own code O(4) change over 0-1000ppm
%         eta = etain.*1e2.*exp(-0.00921.*Xm).*exp(E./(Rg.*T)); %8.73e16
%         
%     case 16 % Our own code O(10) change over 0-1000ppm
%         eta = etain.*1e5.*exp(-0.02303.*Xm).*exp(E./(Rg.*T)); %8.73e19
%         
%     case 17 % Our own code O(0) change over 0-1000ppm
%         eta = etain.*exp(0.*Xm).*exp(E./(Rg.*T)); %8.73e14
%         
%     case 20
%         eta = 8.73e14.*exp(E./(Rg.*T));
%         
%     case 21
%         eta = etain.*exp(-0.004605.*Xm).*exp(E./(Rg.*T));