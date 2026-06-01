%% % add to analytical solution
% sup figures: temporal evolution of piecewise solutions

%% data
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

R0Pu_Ur=0.0068;
lamU=Xe_input.lamUr; % [yr-1]
lamPu=Xe_input.lamPu; % [yr-1]
tau_U=1/lamU/1e9; % [yr]->[Ga]
tau_Pu=1/lamPu/1e9; % [yr]->[Ga]

%% input

tpd=4.6; % Gyr

% Xepd=5e5; % atoms/g; 4.3e5-9.2e5
Xepd=4.3e5;
Xei=3.2e8; % 3.2e7-3.2e8
alpha = Xepd/Xei;

ri=Xeren; % 128Xe/130Xe of initial mantle
ra=Xerpa; % 128Xe/130Xe of present day atmosphere
rpd=Xe_obs.Xe128v130(2);

rfi=Xe132ren;
rfa=Xe132rpa;
rfa1=Xe132rpa;
rfa2=Xe132rpa;
rfpd=1/Xe_obs.Xe130v132(1);
YPu=YPu132;
YU=YUr132;


%%
k128=(ri-rpd)/(ri-ra);
tauD0=tpd/log(Xei/ (Xepd.* (1-k128) ) );
KPu0=YPu./(tau_Pu./tauD0-1);
KU0=YU./(tau_U./tauD0-1);
UAD0=KPu0 .* R0Pu_Ur .* ( exp(-tpd./tau_Pu) - exp(-tpd./tauD0) ) ...
     +KU0 .* ( exp(-tpd./tau_U) - exp(-tpd./tauD0) );
UBSE0=Xepd.* ...
    ( rfpd - rfa.*k128 -rfi.* (1-k128) ) ...
    ./UAD0 ...
    .*exp(-tpd./tau_U)...
    .*238/NA*1e9;


%%

%
n=0:0.01:1;
beta=1; % at T, how much fraction of U left in mantle after CC extraction

%
rT=( (1-n)*ri + n.*(1-alpha)*ra ) ./ (1-n.*alpha);
x1=(1/alpha-n)./(1-n);
x2=(rT-ra)./(rpd-ra);

%
TD=n*tpd.*log(x1)./log(x1.^n.*x2);
tauDD=tpd./log(x1.^n.*x2);
tauDD1=n.*tauDD;
tauDD2=tauDD;

KPu1=YPu./(tau_Pu./tauDD1-1);
KU1=YU./(tau_U./tauDD1-1);
KPu2=YPu./(tau_Pu./tauDD2-1);
KU2=YU./(tau_U./tauDD2-1);

RTPu_Ur = R0Pu_Ur .* exp( -TD./tau_Pu + TD./tau_U );

UAD1=KPu1 .* R0Pu_Ur .* ( exp(-TD./tau_Pu) - exp(-TD./tauDD1) ) ...
     +KU1 .* ( exp(-TD./tau_U) - exp(-TD./tauDD1) );
UAD2=KPu2 .* RTPu_Ur .* ( exp(-(tpd-TD)./tau_Pu) - exp(-(tpd-TD)./tauDD2) ) ...
     +KU2 .* ( exp(-(tpd-TD)./tau_U) - exp(-(tpd-TD)./tauDD2) );

rfTD = ( ...
    UAD1.* exp(TD/tau_U)./ beta.* ( rfpd - rfa2.*( 1 - exp( -(tpd-TD)./tauDD2 ) ) ) ...
    +UAD2.* ( n.*rfa1.* ( 1 - exp( -TD./tauDD1 ) ) + rfi./alpha .* exp( -TD./tauDD1 )  )...
    )./ ( UAD2+ UAD1.* exp( TD./tau_U- (tpd-TD)./tauDD2 )./ beta );


UBSED=Xepd.* ...
    ( rfTD - n.*rfa1.*( 1- exp(-TD./tauDD1 ) ) -rfi./alpha .* exp( -TD./tauDD1 ) ) ...
    ./UAD1 ...
    .*exp(-tpd./tau_U)...
    .*238/NA*1e9;


%
%
TR=tpd*log(x1)./log(x1.*x2);
tauDR=tpd./log(x1.*x2);
tauDR1=tauDR;
tauDR2=tauDR;

KPu1=YPu./(tau_Pu./tauDR1-1);
KU1=YU./(tau_U./tauDR1-1);
KPu2=YPu./(tau_Pu./tauDR2-1);
KU2=YU./(tau_U./tauDR2-1);

RTPu_Ur = R0Pu_Ur .* exp( -TR./tau_Pu + TR./tau_U );

UAR1=KPu1 .* R0Pu_Ur .* ( exp(-TR./tau_Pu) - exp(-TR./tauDR1) ) ...
     +KU1 .* ( exp(-TR./tau_U) - exp(-TR./tauDR1) );
UAR2=KPu2 .* RTPu_Ur .* ( exp(-(tpd-TR)./tau_Pu) - exp(-(tpd-TR)./tauDR2) ) ...
     +KU2 .* ( exp(-(tpd-TR)./tau_U) - exp(-(tpd-TR)./tauDR2) );

rfTR = ( ...
    UAR1.* exp(TR/tau_U)./ beta .* ( rfpd - rfa2.*( 1 - exp( -(tpd-TR)./tauDR2 ) ) ) ...
    +UAR2.* ( n.*rfa1.* ( 1 - exp( -TR./tauDR1 ) ) + rfi./alpha .* exp( -TR./tauDR1 )  )...
    )./ ( UAR2+ UAR1.* exp( TR./tau_U- (tpd-TR)./tauDR2 )./ beta );


UBSER=Xepd.* ...
    ( rfTR - n.*rfa1.*( 1- exp(-TR./tauDR1 ) ) -rfi./alpha .* exp( -TR./tauDR1 ) ) ...
    ./UAR1 ...
    .*exp(-tpd./tau_U)...
    .*238/NA*1e9;



%%
%
YPu130=0;YUr130=0;
Xe130_i=Xei;

t=0:0.01:4.6;

%% constant D&R solution
Xeslab=8e6;
D=M/tauD0; % 5.128
R=k128*Xepd/Xeslab*D;
%
% U=
U238bse=UBSE0*1e-9; % [ppb]->[g/g]
U_i=U238bse*exp(Xe_input.lamUr*t_pd/yr_s)/238*NA; % [g/g]->[atmos/g]
Pu_i=U_i*R0Pu_Ur; % [atmos/g]

% calculate dependent params
Xe130_R=Xeslab*R/D;
% pd atm
Xe128_R1=Xe130_R*Xerpa;
Xe132_R1=Xe130_R*Xe132rpa;

tau_D=M/D; % [Ga]

KPu130=0;
KU130=0;
KPu128=0;
KU128=0;
KPu132=YPu132/(tau_Pu/tau_D-1);
KU132=YUr132/(tau_U/tau_D-1);

Xe_R=Xe130_R;Xe_i=Xe130_i;
KPu=KPu130;KU=KU130;
Xe130=(Xe_i-Xe_R-KPu*Pu_i-KU*U_i)*exp(-t./tau_D)+Xe_R...
        +KPu*Pu_i*exp(-t./tau_Pu)+KU*U_i*exp(-t./tau_U);

Xe_R=Xe128_R1;Xe_i=Xe130_i*ri;
KPu=KPu128;KU=KU128;
Xe128=(Xe_i-Xe_R-KPu*Pu_i-KU*U_i)*exp(-t./tau_D)+Xe_R...
        +KPu*Pu_i*exp(-t./tau_Pu)+KU*U_i*exp(-t./tau_U);

Xe_R=Xe132_R1;Xe_i=Xe130_i*rfi;
KPu=KPu132;KU=KU132;
Xe132=(Xe_i-Xe_R-KPu*Pu_i-KU*U_i)*exp(-t./tau_D)+Xe_R...
        +KPu*Pu_i*exp(-t./tau_Pu)+KU*U_i*exp(-t./tau_U);


%% step D solution

kd=35; % n=0.34
nd=n(kd);Xe_R1=nd*Xepd;Xe_R2=Xepd;Xe_T=Xepd;
tau_D1=tauDD1(kd);tau_D2=tauDD2(kd);TT=TD(kd);
U_i=UBSED(kd)/1e9/238*NA.*exp(tpd./tau_U);Pu_i=U_i*R0Pu_Ur;
U_T=U_i.*exp(-TT./tau_U)*beta;Pu_T=Pu_i.*exp(-TT./tau_Pu)*beta;


Xe_i=Xe130_i;YPu=YPu130;YU=YUr130;

KPu1=YPu./(tau_Pu./tau_D1-1);
KU1=YU./(tau_U./tau_D1-1);
KPu2=YPu./(tau_Pu./tau_D2-1);
KU2=YU./(tau_U./tau_D2-1);

Xe130D=(t<TT).*...
    ( (Xe_i-Xe_R1-KPu1*Pu_i-KU1*U_i).*exp(-t./tau_D1)+Xe_R1...
        +KPu1*Pu_i*exp(-t./tau_Pu)+KU1*U_i*exp(-t./tau_U) )...
        +(t>TT).*...
        ( (Xe_T-Xe_R2-KPu2*Pu_T-KU2*U_T)*exp(-(t-TT)./tau_D2)+Xe_R2...
        +KPu2*Pu_T*exp(-(t-TT)./tau_Pu)+KU2*U_T*exp(-(t-TT)./tau_U) );

%

t=0:0.01:4.6;

kr=1;
nr=n(kr);Xe_R1=nr*Xepd;Xe_R2=Xepd;Xe_T=Xepd;
tau_D1=tauDR1(kr);tau_D2=tauDR2(kr);TT=TR(kr);
U_i=UBSER(kr)/1e9/238*NA.*exp(tpd./tau_U);Pu_i=U_i*R0Pu_Ur;
U_T=U_i.*exp(-TT./tau_U)*beta;Pu_T=Pu_i.*exp(-TT./tau_Pu)*beta;


Xe_i=Xe130_i;YPu=YPu130;YU=YUr130;

KPu1=YPu./(tau_Pu./tau_D1-1);
KU1=YU./(tau_U./tau_D1-1);
KPu2=YPu./(tau_Pu./tau_D2-1);
KU2=YU./(tau_U./tau_D2-1);

Xe130R=(t<TT).*...
    ( (Xe_i-Xe_R1-KPu1*Pu_i-KU1*U_i).*exp(-t./tau_D1)+Xe_R1...
        +KPu1*Pu_i*exp(-t./tau_Pu)+KU1*U_i*exp(-t./tau_U) )...
        +(t>TT).*...
        ( (Xe_T-Xe_R2-KPu2*Pu_T-KU2*U_T)*exp(-(t-TT)./tau_D2)+Xe_R2...
        +KPu2*Pu_T*exp(-(t-TT)./tau_Pu)+KU2*U_T*exp(-(t-TT)./tau_U) );

plot(t,Xe130D);hold on;plot(t,Xe130R);

%%

% Xe_R=Xe128_R1;Xe_i=Xe128_i;
% KPu=KPu128;KU=KU128;
% Xe128=(Xe_i-Xe_R-KPu*Pu_i-KU*U_i)*exp(-t./tau_D)+Xe_R...
%         +KPu*Pu_i*exp(-t./tau_Pu)+KU*U_i*exp(-t./tau_U);

%
YPu128=0;YUr128=0;
Xe128_i=Xei*ri;

Xe_i=Xe128_i;YPu=YPu128;YU=YUr128;
t=0:0.01:4.6;

kd=35;
nd=n(kd);rTT=rT(kd);
Xe_R1=nd*Xepd*ra;Xe_R2=Xepd*ra;Xe_T=Xepd*rTT;
tau_D1=tauDD1(kd);tau_D2=tauDD2(kd);TT=TD(kd);
U_i=UBSED(kd)/1e9/238*NA.*exp(tpd./tau_U);Pu_i=U_i*R0Pu_Ur;
U_T=U_i.*exp(-TT./tau_U)*beta;Pu_T=Pu_i.*exp(-TT./tau_Pu)*beta;

KPu1=YPu./(tau_Pu./tau_D1-1);
KU1=YU./(tau_U./tau_D1-1);
KPu2=YPu./(tau_Pu./tau_D2-1);
KU2=YU./(tau_U./tau_D2-1);

Xe128D=(t<TT).*...
    ( (Xe_i-Xe_R1-KPu1*Pu_i-KU1*U_i).*exp(-t./tau_D1)+Xe_R1...
        +KPu1*Pu_i*exp(-t./tau_Pu)+KU1*U_i*exp(-t./tau_U) )...
        +(t>TT).*...
        ( (Xe_T-Xe_R2-KPu2*Pu_T-KU2*U_T)*exp(-(t-TT)./tau_D2)+Xe_R2...
        +KPu2*Pu_T*exp(-(t-TT)./tau_Pu)+KU2*U_T*exp(-(t-TT)./tau_U) );

%

kr=1;
nr=n(kr);rTT=rT(kr);
Xe_R1=nr*Xepd*ra;Xe_R2=Xepd*ra;Xe_T=Xepd*rTT;
tau_D1=tauDR1(kr);tau_D2=tauDR2(kr);TT=TR(kr);
U_i=UBSER(kr)/1e9/238*NA.*exp(tpd./tau_U);Pu_i=U_i*R0Pu_Ur;
U_T=U_i.*exp(-TT./tau_U)*beta;Pu_T=Pu_i.*exp(-TT./tau_Pu)*beta;


KPu1=YPu./(tau_Pu./tau_D1-1);
KU1=YU./(tau_U./tau_D1-1);
KPu2=YPu./(tau_Pu./tau_D2-1);
KU2=YU./(tau_U./tau_D2-1);

Xe128R=(t<TT).*...
    ( (Xe_i-Xe_R1-KPu1*Pu_i-KU1*U_i).*exp(-t./tau_D1)+Xe_R1...
        +KPu1*Pu_i*exp(-t./tau_Pu)+KU1*U_i*exp(-t./tau_U) )...
        +(t>TT).*...
        ( (Xe_T-Xe_R2-KPu2*Pu_T-KU2*U_T)*exp(-(t-TT)./tau_D2)+Xe_R2...
        +KPu2*Pu_T*exp(-(t-TT)./tau_Pu)+KU2*U_T*exp(-(t-TT)./tau_U) );

plot(t,Xe128D./Xe130D);hold on;plot(t,Xe128R./Xe130R);


%%
% Xe_R=Xe132_R1;Xe_i=Xe132_i;
% KPu=KPu132;KU=KU132;
% Xe132=(Xe_i-Xe_R-KPu*Pu_i-KU*U_i)*exp(-t./tau_D)+Xe_R...
%         +KPu*Pu_i*exp(-t./tau_Pu)+KU*U_i*exp(-t./tau_U);


Xe132_i=Xei*rfi;
Xe_i=Xe132_i;YPu=YPu132;YU=YUr132;
t=0:0.01:4.6;

kd=35;
nd=n(kd);rTT=rfTD(kd);
Xe_R1=nd*Xepd*rfa1;Xe_R2=Xepd*rfa2;Xe_T=Xepd*rTT;
tau_D1=tauDD1(kd);tau_D2=tauDD2(kd);TT=TD(kd);
U_i=UBSED(kd)/1e9/238*NA.*exp(tpd./tau_U);Pu_i=U_i*R0Pu_Ur;
U_T=U_i.*exp(-TT./tau_U)*beta;Pu_T=Pu_i.*exp(-TT./tau_Pu)*beta;

KPu1=YPu./(tau_Pu./tau_D1-1);
KU1=YU./(tau_U./tau_D1-1);
KPu2=YPu./(tau_Pu./tau_D2-1);
KU2=YU./(tau_U./tau_D2-1);

Xe132D=(t<TT).*...
    ( (Xe_i-Xe_R1-KPu1*Pu_i-KU1*U_i).*exp(-t./tau_D1)+Xe_R1...
        +KPu1*Pu_i*exp(-t./tau_Pu)+KU1*U_i*exp(-t./tau_U) )...
        +(t>TT).*...
        ( (Xe_T-Xe_R2-KPu2*Pu_T-KU2*U_T)*exp(-(t-TT)./tau_D2)+Xe_R2...
        +KPu2*Pu_T*exp(-(t-TT)./tau_Pu)+KU2*U_T*exp(-(t-TT)./tau_U) );

U_i=U238bse*exp(Xe_input.lamUr*t_pd/yr_s)/238*NA;Pu_i=U_i*R0Pu_Ur;
U_T=U_i.*exp(-TT./tau_U)*beta;Pu_T=Pu_i.*exp(-TT./tau_Pu)*beta;
Xe132D_U=(t<TT).*...
    ( (Xe_i-Xe_R1-KPu1*Pu_i-KU1*U_i).*exp(-t./tau_D1)+Xe_R1...
        +KPu1*Pu_i*exp(-t./tau_Pu)+KU1*U_i*exp(-t./tau_U) )...
        +(t>TT).*...
        ( (Xe_T-Xe_R2-KPu2*Pu_T-KU2*U_T)*exp(-(t-TT)./tau_D2)+Xe_R2...
        +KPu2*Pu_T*exp(-(t-TT)./tau_Pu)+KU2*U_T*exp(-(t-TT)./tau_U) );

%

kr=1;
nr=n(kr);rTT=rfTR(kr);
Xe_R1=nr*Xepd*rfa1;Xe_R2=Xepd*rfa2;Xe_T=Xepd*rTT;
tau_D1=tauDR1(kr);tau_D2=tauDR2(kr);TT=TR(kr);
U_i=UBSER(kr)/1e9/238*NA.*exp(tpd./tau_U);Pu_i=U_i*R0Pu_Ur;
U_T=U_i.*exp(-TT./tau_U)*beta;Pu_T=Pu_i.*exp(-TT./tau_Pu)*beta;


KPu1=YPu./(tau_Pu./tau_D1-1);
KU1=YU./(tau_U./tau_D1-1);
KPu2=YPu./(tau_Pu./tau_D2-1);
KU2=YU./(tau_U./tau_D2-1);

Xe132R=(t<TT).*...
    ( (Xe_i-Xe_R1-KPu1*Pu_i-KU1*U_i).*exp(-t./tau_D1)+Xe_R1...
        +KPu1*Pu_i*exp(-t./tau_Pu)+KU1*U_i*exp(-t./tau_U) )...
        +(t>TT).*...
        ( (Xe_T-Xe_R2-KPu2*Pu_T-KU2*U_T)*exp(-(t-TT)./tau_D2)+Xe_R2...
        +KPu2*Pu_T*exp(-(t-TT)./tau_Pu)+KU2*U_T*exp(-(t-TT)./tau_U) );

plot(t,Xe132D./Xe130D);hold on;plot(t,Xe132R./Xe130R);

%%
% Xeslab=8e6;

kd=35;nd=n(kd);
%Xe_R1=nd*Xepd;Xe_R2=Xepd;
%Xe_T=Xepd;
tau_D1=tauDD1(kd);tau_D2=tauDD2(kd);TT=TD(kd);
% tauDplot=(t<TT).*tau_D1+(t>=TT).*tau_D2;
D=M/tauD0;
D1=M/tau_D1;D2=M/tau_D2;
R2=Xepd/Xeslab*D2;R1=R2;
DDplot=(t<TT).*D1+(t>=TT).*D2;
RDplot=(t<TT).*R1+(t>=TT).*R2;

D0plot=(t<TT).*D+(t>=TT).*D;
R0plot=(t<TT).*R+(t>=TT).*R;

kr=1;nr=n(kr);

Xe_R1=nr*Xepd;Xe_R2=Xepd;
% Xe_T=Xepd*rTT;
% tau_D1=tauDR1(kr);tau_D2=tauDR2(kr);
TT=TR(kr);
tau_D00=tauDR(kr);
D1=M/tau_D00;D2=M/tau_D00;
R2=Xepd/Xeslab*D2;R1=nr*R2;
R=k128*Xepd/Xeslab*D;
DRplot=(t<TT).*D1+(t>=TT).*D2;
RRplot=(t<TT).*R1+(t>=TT).*R2;


%%
close all;
figure;
set(gcf,'color','white','Units','centimeters','Position',[2 5 17.8 12]);

gca=subplot(3,2,1);
% semilogy(t,DDplot,'LineWidth',2);hold on;
% semilogy(t,DRplot,'LineWidth',2);
plot(t,DDplot,'LineWidth',2,'color','r');hold on;
plot(t,DRplot,'LineWidth',2,'color','b');hold on;
plot(t,D0plot,'LineWidth',2,'color','k');
xlabel('Time (Gyr)');
ylabel('Degessing flux (kg/Gyr)');
set(gca,'FontUnits','points','FontSize',7,'FontWeight','bold',...
    'LineWidth',1.5,'LabelFontSizeMultiplier',1);
legend('n=0.34, m=1',...
    'n=1, m=0','n=1, m=1',...
    'location',[0.177,0.179,0.19,0.0853],'fontsize',10);
text(-0.126,1.01,'a','Units','normalized','FontWeight','bold','FontSize',8);

gca=subplot(3,2,3);
plot(t,RDplot,'LineWidth',2,'color','r');hold on;
plot(t,RRplot,'LineWidth',2,'color','b');hold on;
plot(t,R0plot,'LineWidth',2,'color','k');
% legend('piecewise D',...
%     'piecewise R');
xlabel('Time (Gyr)');
ylabel('Regessing flux (kg/Gyr)');
set(gca,'FontUnits','points','FontSize',7,'FontWeight','bold',...
    'LineWidth',1.5,'LabelFontSizeMultiplier',1);
xlim([0 4.6]);
text(-0.126,1.01,'b','Units','normalized','FontWeight','bold','FontSize',8);


%

gca=subplot(3,2,2);
plot(t,Xe130D,'LineWidth',2,'color','r');hold on;
plot(t,Xe130R,'LineWidth',2,'color','b');hold on;
plot(t,Xe130,'LineWidth',2,'color','k');
% legend('piecewise D','piecewise R','constant D&R');
xlabel('Time (Gyr)');
ylabel('^{130}Xe');
set(gca,'FontUnits','points','FontSize',7,'FontWeight','bold',...
    'LineWidth',1.5,'LabelFontSizeMultiplier',1);
xlim([0 4.6]);
text(-0.126,1.01,'c','Units','normalized','FontWeight','bold','FontSize',8);

gca=subplot(3,2,4);
plot(t,Xe128D./Xe130D,'LineWidth',2,'Color','r');hold on;
plot(t,Xe128R./Xe130R,'LineWidth',2,'Color','b');hold on;
plot(t,Xe128./Xe130,'LineWidth',2,'Color','k');
% legend('piecewise D',...
%     'piecewise R');
xlabel('Time (Gyr)');
ylabel('^{128}Xe/^{130}Xe');
set(gca,'FontUnits','points','FontSize',7,'FontWeight','bold',...
    'LineWidth',1.5,'LabelFontSizeMultiplier',1);
xlim([0 4.6]);
text(-0.126,1.01,'d','Units','normalized','FontWeight','bold','FontSize',8);


gca=subplot(3,2,6);
plot(t,Xe132D./Xe130D,'LineWidth',2,'color','r');hold on;
plot(t,Xe132R./Xe130R,'LineWidth',2,'color','b');hold on;
plot(t,Xe132./Xe130,'LineWidth',2,'color','k');
% legend('piecewise D',...
%     'piecewise R');
xlabel('Time (Gyr)');
ylabel('^{132}Xe/^{130}Xe');
set(gca,'FontUnits','points','FontSize',7,'FontWeight','bold',...
    'LineWidth',1.5,'LabelFontSizeMultiplier',1);
xlim([0 4.6]);
text(-0.126,1.01,'e','Units','normalized','FontWeight','bold','FontSize',8);

% exportgraphics(gcf, 'test_figS1s1_analytical_step2_3e8.pdf', ...
%     'ContentType','image', ...   % vector/image
%     'BackgroundColor','white', ... % 
%     'Resolution',600); 



