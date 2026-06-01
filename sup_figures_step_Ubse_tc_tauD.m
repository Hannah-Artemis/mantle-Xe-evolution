%% sup_figures_step_new



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
Xe_obs1 = struct( 'Xe',       [4.3e5	9.2e5], ... % atoms/g
                 'Xe128v130',     [0.475	0.478], ...
                 'Xe128v132',     [0.069	0.071], ...
                 'Xe132v130',     [1/0.1493 1/0.1445], ...
                 'Xe131v130',     [0.7608/0.1493	0.7786/0.1445], ...
                 'Xe134v130',     [0.4082/0.1493	0.4302/0.1445], ...
                 'Xe136v130',     [0.3559/0.1493	0.3835/0.1445] ...
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

m128_range_cc=[(Xercc-Xe_obs.Xe128v130(2))/(Xercc-Xerpa),(Xercc-Xe_obs.Xe128v130(1))/(Xercc-Xerpa)];
m128_range_en=[(Xeren-Xe_obs.Xe128v130(2))/(Xeren-Xerpa),(Xeren-Xe_obs.Xe128v130(1))/(Xeren-Xerpa)];


%% set
%
tpd=4.6; % Gyr

Xepdmin=Xe_obs.Xe(1);
Xepdmax=Xe_obs.Xe(2);
Xeimin=3.2e7;
Xeimax=3.2e8;
% k128min=m128_range_en(1);
% k128max=m128_range_en(2);
rpdmin=Xe_obs.Xe128v130(1);
rpdmax=Xe_obs.Xe128v130(2);

%
id=0;
input_test=zeros(1000,3);
for Xepd=Xepdmin:(Xepdmax-Xepdmin)/9:Xepdmax
    for Xei=Xeimin:(Xeimax-Xeimin)/9:Xeimax
        for rpd=rpdmin:(rpdmax-rpdmin)/9:rpdmax
            % for n=0.01:0.0098:0.99
                id=id+1;
                input_test(id,1)=Xepd;
                input_test(id,2)=Xei;
                input_test(id,3)=rpd;
                % input_test(id,4)=n;
            % end      
        end
    end
end

Xepd=input_test(:,1);
Xei=input_test(:,2);
rpd=input_test(:,3);
% n=input_test(:,4);

alpha = Xepd./Xei;

ri=Xeren; % 128Xe/130Xe of initial mantle
ra=Xerpa; % 128Xe/130Xe of present day atmosphere



%% U from constant D & R
U_constant_extrema=zeros(4,9);



% 131
niso=131;

U_constant_extrema(1,1)=niso;
rfi=Xe131ren;
rfa=Xe131rpa;
rfa1=Xe131rpa;
rfa2=Xe131rpa;
YPu=YPu131;
YU=YUr131;

% min
rpd=rpdmax;
rfpd=Xe_obs1.Xe131v130(1); % min fXe/130Xe

k128=(ri-rpd)./(ri-ra);
tauD0=tpd./log(Xei./ (Xepd.* (1-k128) ) );
KPu0=YPu./(tau_Pu./tauD0-1);
KU0=YU./(tau_U./tauD0-1);
UAD0=KPu0 .* R0Pu_Ur .* ( exp(-tpd./tau_Pu) - exp(-tpd./tauD0) ) ...
     +KU0 .* ( exp(-tpd./tau_U) - exp(-tpd./tauD0) );
UBSE0=Xepd.* ...
    ( rfpd - rfa.*k128 -rfi.* (1-k128) ) ...
    ./UAD0 ...
    .*exp(-tpd./tau_U)...
    .*238/NA*1e9;

[U_constant_extrema(1,2),U_constant_extrema(1,3)]=min(UBSE0);
U_constant_extrema(1,6)=min(UBSE0(Xei==3.2e7));
U_constant_extrema(1,8)=min(UBSE0(Xei==3.2e8));

% max
rpd=rpdmin;
rfpd=Xe_obs1.Xe131v130(2); % min fXe/130Xe

k128=(ri-rpd)./(ri-ra);
tauD0=tpd./log(Xei./ (Xepd.* (1-k128) ) );
KPu0=YPu./(tau_Pu./tauD0-1);
KU0=YU./(tau_U./tauD0-1);
UAD0=KPu0 .* R0Pu_Ur .* ( exp(-tpd./tau_Pu) - exp(-tpd./tauD0) ) ...
     +KU0 .* ( exp(-tpd./tau_U) - exp(-tpd./tauD0) );
UBSE0=Xepd.* ...
    ( rfpd - rfa.*k128 -rfi.* (1-k128) ) ...
    ./UAD0 ...
    .*exp(-tpd./tau_U)...
    .*238/NA*1e9;

[U_constant_extrema(1,4),U_constant_extrema(1,5)]=max(UBSE0);
U_constant_extrema(1,7)=max(UBSE0(Xei==3.2e7));
U_constant_extrema(1,9)=max(UBSE0(Xei==3.2e8));


% 132
niso=132;

U_constant_extrema(2,1)=niso;
rfi=Xe132ren;
rfa=Xe132rpa;
rfa1=Xe132rpa;
rfa2=Xe132rpa;
YPu=YPu132;
YU=YUr132;

% min
rpd=rpdmax;
rfpd=Xe_obs1.Xe132v130(1); % min fXe/130Xe

k128=(ri-rpd)./(ri-ra);
tauD0=tpd./log(Xei./ (Xepd.* (1-k128) ) );
KPu0=YPu./(tau_Pu./tauD0-1);
KU0=YU./(tau_U./tauD0-1);
UAD0=KPu0 .* R0Pu_Ur .* ( exp(-tpd./tau_Pu) - exp(-tpd./tauD0) ) ...
     +KU0 .* ( exp(-tpd./tau_U) - exp(-tpd./tauD0) );
UBSE0=Xepd.* ...
    ( rfpd - rfa.*k128 -rfi.* (1-k128) ) ...
    ./UAD0 ...
    .*exp(-tpd./tau_U)...
    .*238/NA*1e9;

[U_constant_extrema(2,2),U_constant_extrema(2,3)]=min(UBSE0);
U_constant_extrema(2,6)=min(UBSE0(Xei==3.2e7));
U_constant_extrema(2,8)=min(UBSE0(Xei==3.2e8));

% max
rpd=rpdmin;
rfpd=Xe_obs1.Xe132v130(2); % min fXe/130Xe

k128=(ri-rpd)./(ri-ra);
tauD0=tpd./log(Xei./ (Xepd.* (1-k128) ) );
KPu0=YPu./(tau_Pu./tauD0-1);
KU0=YU./(tau_U./tauD0-1);
UAD0=KPu0 .* R0Pu_Ur .* ( exp(-tpd./tau_Pu) - exp(-tpd./tauD0) ) ...
     +KU0 .* ( exp(-tpd./tau_U) - exp(-tpd./tauD0) );
UBSE0=Xepd.* ...
    ( rfpd - rfa.*k128 -rfi.* (1-k128) ) ...
    ./UAD0 ...
    .*exp(-tpd./tau_U)...
    .*238/NA*1e9;

[U_constant_extrema(2,4),U_constant_extrema(2,5)]=max(UBSE0);
U_constant_extrema(2,7)=max(UBSE0(Xei==3.2e7));
U_constant_extrema(2,9)=max(UBSE0(Xei==3.2e8));

% 134
niso=134;

U_constant_extrema(3,1)=niso;
rfi=Xe134ren;
rfa=Xe134rpa;
rfa1=Xe134rpa;
rfa2=Xe134rpa;
YPu=YPu134;
YU=YUr134;

% min
rpd=rpdmax;
rfpd=Xe_obs1.Xe134v130(1); % min fXe/130Xe

k128=(ri-rpd)./(ri-ra);
tauD0=tpd./log(Xei./ (Xepd.* (1-k128) ) );
KPu0=YPu./(tau_Pu./tauD0-1);
KU0=YU./(tau_U./tauD0-1);
UAD0=KPu0 .* R0Pu_Ur .* ( exp(-tpd./tau_Pu) - exp(-tpd./tauD0) ) ...
     +KU0 .* ( exp(-tpd./tau_U) - exp(-tpd./tauD0) );
UBSE0=Xepd.* ...
    ( rfpd - rfa.*k128 -rfi.* (1-k128) ) ...
    ./UAD0 ...
    .*exp(-tpd./tau_U)...
    .*238/NA*1e9;

[U_constant_extrema(3,2),U_constant_extrema(3,3)]=min(UBSE0);
U_constant_extrema(3,6)=min(UBSE0(Xei==3.2e7));
U_constant_extrema(3,8)=min(UBSE0(Xei==3.2e8));

% max
rpd=rpdmin;
rfpd=Xe_obs1.Xe134v130(2); % min fXe/130Xe

k128=(ri-rpd)./(ri-ra);
tauD0=tpd./log(Xei./ (Xepd.* (1-k128) ) );
KPu0=YPu./(tau_Pu./tauD0-1);
KU0=YU./(tau_U./tauD0-1);
UAD0=KPu0 .* R0Pu_Ur .* ( exp(-tpd./tau_Pu) - exp(-tpd./tauD0) ) ...
     +KU0 .* ( exp(-tpd./tau_U) - exp(-tpd./tauD0) );
UBSE0=Xepd.* ...
    ( rfpd - rfa.*k128 -rfi.* (1-k128) ) ...
    ./UAD0 ...
    .*exp(-tpd./tau_U)...
    .*238/NA*1e9;

[U_constant_extrema(3,4),U_constant_extrema(3,5)]=max(UBSE0);
U_constant_extrema(3,7)=max(UBSE0(Xei==3.2e7));
U_constant_extrema(3,9)=max(UBSE0(Xei==3.2e8));

% 136
niso=136;

U_constant_extrema(4,1)=niso;
rfi=Xe136ren;
rfa=Xe136rpa;
rfa1=Xe136rpa;
rfa2=Xe136rpa;
YPu=YPu;
YU=YUr;

% min
rpd=rpdmax;
rfpd=Xe_obs1.Xe136v130(1); % min fXe/130Xe

k128=(ri-rpd)./(ri-ra);
tauD0=tpd./log(Xei./ (Xepd.* (1-k128) ) );
KPu0=YPu./(tau_Pu./tauD0-1);
KU0=YU./(tau_U./tauD0-1);
UAD0=KPu0 .* R0Pu_Ur .* ( exp(-tpd./tau_Pu) - exp(-tpd./tauD0) ) ...
     +KU0 .* ( exp(-tpd./tau_U) - exp(-tpd./tauD0) );
UBSE0=Xepd.* ...
    ( rfpd - rfa.*k128 -rfi.* (1-k128) ) ...
    ./UAD0 ...
    .*exp(-tpd./tau_U)...
    .*238/NA*1e9;

[U_constant_extrema(4,2),U_constant_extrema(4,3)]=min(UBSE0);
U_constant_extrema(4,6)=min(UBSE0(Xei==3.2e7));
U_constant_extrema(4,8)=min(UBSE0(Xei==3.2e8));

% max
rpd=rpdmin;
rfpd=Xe_obs1.Xe136v130(2); % min fXe/130Xe

k128=(ri-rpd)./(ri-ra);
tauD0=tpd./log(Xei./ (Xepd.* (1-k128) ) );
KPu0=YPu./(tau_Pu./tauD0-1);
KU0=YU./(tau_U./tauD0-1);
UAD0=KPu0 .* R0Pu_Ur .* ( exp(-tpd./tau_Pu) - exp(-tpd./tauD0) ) ...
     +KU0 .* ( exp(-tpd./tau_U) - exp(-tpd./tauD0) );
UBSE0=Xepd.* ...
    ( rfpd - rfa.*k128 -rfi.* (1-k128) ) ...
    ./UAD0 ...
    .*exp(-tpd./tau_U)...
    .*238/NA*1e9;

[U_constant_extrema(4,4),U_constant_extrema(4,5)]=max(UBSE0);
U_constant_extrema(4,7)=max(UBSE0(Xei==3.2e7));
U_constant_extrema(4,9)=max(UBSE0(Xei==3.2e8));

UBSE0_m=zeros(6,1);
idx_iso=[131 132 134 136];
UBSE0_m(1)=max(U_constant_extrema(:,2));
UBSE0_m(2)=min(U_constant_extrema(:,4));
UBSE0_m(3)=max(U_constant_extrema(:,6));
UBSE0_m(4)=min(U_constant_extrema(:,7));
UBSE0_m(5)=max(U_constant_extrema(:,8));
UBSE0_m(6)=min(U_constant_extrema(:,9));


%% update set for changing D&R
%
id=0;
input_test_step=zeros(11907,5); %3^3*21*21 
for Xepd=Xepdmin:(Xepdmax-Xepdmin)/2:Xepdmax
    for Xei=Xeimin:(Xeimax-Xeimin)/2:Xeimax
        for rpd=rpdmin:(rpdmax-rpdmin)/2:rpdmax
            %for n=0.01:0.0098:0.99
            for n=0.01:0.049:0.99
                for m=0.01:0.049:0.99
                    id=id+1;
                    input_test_step(id,1)=Xepd;
                    input_test_step(id,2)=Xei;
                    input_test_step(id,3)=rpd;
                    input_test_step(id,4)=n;
                    input_test_step(id,5)=m;
                end
            end      
        end
    end
end

Xepd=input_test_step(:,1);
Xei=input_test_step(:,2);
rpd=input_test_step(:,3);
n=input_test_step(:,4);
m=input_test_step(:,5);


%% Transition time and Degassing timescale with R1=mR2; D2=nD1

%
ri=Xeren; % 128Xe/130Xe of initial mantle
ra=Xerpa; % 128Xe/130Xe of present day atmosphere

Xepd=input_test_step(7,1);
Xei=input_test_step(7,2);
rpd=input_test_step(7,3);

alpha=Xepd./Xei;
k128=(ri-rpd)./(ri-ra);
rT=( (1-n.*m).*ri + n.*m.*(1-alpha).*ra ) ./ (1-n.*m.*alpha);
x1=(1./alpha-n.*m)./(1-n.*m);
x2=(rT-ra)./(rpd-ra);

%
TD=n.*tpd.*log(x1)./log(x1.^n.*x2);
tauDD=tpd./log(x1.^n.*x2);
% tauDD_test=tpd./log(x1.^(n-1)./alpha./(1-k128))
tauDD1=n.*tauDD;
tauDD2=tauDD;

%% collect  extrema T & tau_D data for piecewise D&R with given n&m
%
Xepd=input_test_step(:,1);
Xei=input_test_step(:,2);
rpd=input_test_step(:,3);
n=input_test_step(:,4);
m=input_test_step(:,5);
alpha=Xepd./Xei;


%
Nn=length(0.01:0.049:0.99);
Nm=length(0.01:0.049:0.99);
nset=m(1:Nm);
mset=m(1:Nm);
% Nniso=4; % number of isotope

%
tauDD1_m=zeros(Nn,Nm,10);
tauDD2_m=zeros(Nn,Nm,10);
TD_m=zeros(Nn,Nm,10);

%
for kn=1:Nn
    ntest=nset(kn);
    

    for km=1:Nm
        mtest=mset(km);
        idx_all = find(n == ntest & m==mtest);
        %
        tauDD1_m(kn,km,1)=ntest;tauDD1_m(kn,km,2)=mtest;
        [tauDD1_m(kn,km,3),idx_sub]=min(tauDD1(idx_all)); % min
        index=idx_all(idx_sub);
        tauDD1_m(kn,km,4)=Xepd(index);tauDD1_m(kn,km,5)=Xei(index);tauDD1_m(kn,km,6)=rpd(index);
        [tauDD1_m(kn,km,7),idx_sub]=max(tauDD1(idx_all)); % max
        index=idx_all(idx_sub);
        tauDD1_m(kn,km,8)=Xepd(index);tauDD1_m(kn,km,9)=Xei(index);tauDD1_m(kn,km,10)=rpd(index);
    
        tauDD2_m(kn,km,1)=ntest;tauDD2_m(kn,km,2)=mtest;
        [tauDD2_m(kn,km,3),idx_sub]=min(tauDD2(idx_all));
        index=idx_all(idx_sub);
        tauDD2_m(kn,km,4)=Xepd(index);tauDD2_m(kn,km,5)=Xei(index);tauDD2_m(kn,km,6)=rpd(index);
        [tauDD2_m(kn,km,7),idx_sub]=max(tauDD2(idx_all));
        index=idx_all(idx_sub);
        tauDD2_m(kn,km,8)=Xepd(index);tauDD2_m(kn,km,9)=Xei(index);tauDD2_m(kn,km,10)=rpd(index);
    
        TD_m(kn,km,1)=ntest;TD_m(kn,km,2)=mtest;
        [TD_m(kn,km,3),idx_sub]=min(TD(idx_all));
        index=idx_all(idx_sub);
        TD_m(kn,km,4)=Xepd(index);TD_m(kn,km,5)=Xei(index);TD_m(kn,km,6)=rpd(index);
        [TD_m(kn,km,7),idx_sub]=max(TD(idx_all));
        index=idx_all(idx_sub);
        TD_m(kn,km,8)=Xepd(index);TD_m(kn,km,9)=Xei(index);TD_m(kn,km,10)=rpd(index);
    
    end
end


%% U by piecewise D (piecewise R does not affect U)
beta=1; % at T, how much fraction of U left in mantle after CC extraction
U_D_nocc_extrema=zeros(length(0.01:0.049:0.99),56);

%
alpha=Xepd./Xei;
rT=( (1-n).*ri + n.*(1-alpha).*ra ) ./ (1-n.*alpha);
x1=(1./alpha-n)./(1-n);
x2=(rT-ra)./(rpd-ra);

%
TD=n.*tpd.*log(x1)./log(x1.^n.*x2);
tauDD=tpd./log(x1.^n.*x2);
tauDD1=n.*tauDD;
tauDD2=tauDD;
RTPu_Ur = R0Pu_Ur .* exp( -TD./tau_Pu + TD./tau_U );

% 131
niso=131;
% U_D_nocc_extrema(1,1)=niso;
rfi=Xe131ren;
rfa=Xe131rpa;
rfa1=Xe131rpa;
rfa2=Xe131rpa;
YPu=YPu131;
YU=YUr131;

KPu1=YPu./(tau_Pu./tauDD1-1);
KU1=YU./(tau_U./tauDD1-1);
KPu2=YPu./(tau_Pu./tauDD2-1);
KU2=YU./(tau_U./tauDD2-1);

UAD1=KPu1 .* R0Pu_Ur .* ( exp(-TD./tau_Pu) - exp(-TD./tauDD1) ) ...
     +KU1 .* ( exp(-TD./tau_U) - exp(-TD./tauDD1) );
UAD2=KPu2 .* RTPu_Ur .* ( exp(-(tpd-TD)./tau_Pu) - exp(-(tpd-TD)./tauDD2) ) ...
     +KU2 .* ( exp(-(tpd-TD)./tau_U) - exp(-(tpd-TD)./tauDD2) );

% min
rfpd=Xe_obs1.Xe131v130(1); % min fXe/130Xe

rfTD = ( ...
    UAD1.* exp(TD./tau_U)./ beta.* ( rfpd - rfa2.*( 1 - exp( -(tpd-TD)./tauDD2 ) ) ) ...
    +UAD2.* ( n.*rfa1.* ( 1 - exp( -TD./tauDD1 ) ) + rfi./alpha .* exp( -TD./tauDD1 )  )...
    )./ ( UAD2+ UAD1.* exp( TD./tau_U- (tpd-TD)./tauDD2 )./ beta );

UBSED=Xepd.* ...
    ( rfTD - n.*rfa1.*( 1- exp(-TD./tauDD1 ) ) -rfi./alpha .* exp( -TD./tauDD1 ) ) ...
    ./UAD1 ...
    .*exp(-tpd./tau_U)...
    .*238/NA*1e9;

% U_D_nocc_extrema(1,2)=min(UBSED);
% U_D_nocc_extrema(1,4)=min(UBSED(Xei==3.2e7));
% U_D_nocc_extrema(1,6)=min(UBSED(Xei==3.2e8));
% U_D_nocc_extrema(1,8)=min(UBSED(n==n(3)));

Nn=length(0.01:0.049:0.99);
nset=m(1:Nn);
for kn=1:Nn
    ntest=nset(kn);
    idx_all = find(n == ntest);
    %
    U_D_nocc_extrema(kn,1)=ntest;
    U_D_nocc_extrema(kn,2)=niso;    
    [U_D_nocc_extrema(kn,3),idx_sub]=min(UBSED(idx_all));
    index=idx_all(idx_sub);
    U_D_nocc_extrema(kn,4)=Xepd(index);U_D_nocc_extrema(kn,5)=Xei(index);U_D_nocc_extrema(kn,6)=rpd(index);
    U_D_nocc_extrema(kn,11)=min(UBSED(n==ntest & Xei==3.2e7));
    U_D_nocc_extrema(kn,13)=min(UBSED(n==ntest & Xei==3.2e8));
end

% max
rfpd=Xe_obs1.Xe131v130(2); % min fXe/130Xe

rfTD = ( ...
    UAD1.* exp(TD./tau_U)./ beta.* ( rfpd - rfa2.*( 1 - exp( -(tpd-TD)./tauDD2 ) ) ) ...
    +UAD2.* ( n.*rfa1.* ( 1 - exp( -TD./tauDD1 ) ) + rfi./alpha .* exp( -TD./tauDD1 )  )...
    )./ ( UAD2+ UAD1.* exp( TD./tau_U- (tpd-TD)./tauDD2 )./ beta );

UBSED=Xepd.* ...
    ( rfTD - n.*rfa1.*( 1- exp(-TD./tauDD1 ) ) -rfi./alpha .* exp( -TD./tauDD1 ) ) ...
    ./UAD1 ...
    .*exp(-tpd./tau_U)...
    .*238/NA*1e9;

% U_D_nocc_extrema(1,3)=max(UBSED);
% U_D_nocc_extrema(1,5)=max(UBSED(Xei==3.2e7));
% U_D_nocc_extrema(1,7)=max(UBSED(Xei==3.2e8));
% U_D_nocc_extrema(1,9)=max(UBSED(n==n(3)));

for kn=1:length(nset)
    ntest=nset(kn);
    idx_all = find(n == ntest);
    %   
    [U_D_nocc_extrema(kn,7),idx_sub]=max(UBSED(idx_all));
    index=idx_all(idx_sub);
    U_D_nocc_extrema(kn,8)=Xepd(index);U_D_nocc_extrema(kn,9)=Xei(index);U_D_nocc_extrema(kn,10)=rpd(index);
    U_D_nocc_extrema(kn,12)=max(UBSED(n==ntest & Xei==3.2e7));
    U_D_nocc_extrema(kn,14)=max(UBSED(n==ntest & Xei==3.2e8));
end

% 132
niso=132;
% U_D_nocc_extrema(2,1)=niso;
rfi=Xe132ren;
rfa=Xe132rpa;
rfa1=Xe132rpa;
rfa2=Xe132rpa;
YPu=YPu132;
YU=YUr132;


KPu1=YPu./(tau_Pu./tauDD1-1);
KU1=YU./(tau_U./tauDD1-1);
KPu2=YPu./(tau_Pu./tauDD2-1);
KU2=YU./(tau_U./tauDD2-1);


UAD1=KPu1 .* R0Pu_Ur .* ( exp(-TD./tau_Pu) - exp(-TD./tauDD1) ) ...
     +KU1 .* ( exp(-TD./tau_U) - exp(-TD./tauDD1) );
UAD2=KPu2 .* RTPu_Ur .* ( exp(-(tpd-TD)./tau_Pu) - exp(-(tpd-TD)./tauDD2) ) ...
     +KU2 .* ( exp(-(tpd-TD)./tau_U) - exp(-(tpd-TD)./tauDD2) );

% min
rfpd=Xe_obs1.Xe132v130(1); % min fXe/130Xe

rfTD = ( ...
    UAD1.* exp(TD./tau_U)./ beta.* ( rfpd - rfa2.*( 1 - exp( -(tpd-TD)./tauDD2 ) ) ) ...
    +UAD2.* ( n.*rfa1.* ( 1 - exp( -TD./tauDD1 ) ) + rfi./alpha .* exp( -TD./tauDD1 )  )...
    )./ ( UAD2+ UAD1.* exp( TD./tau_U- (tpd-TD)./tauDD2 )./ beta );

UBSED=Xepd.* ...
    ( rfTD - n.*rfa1.*( 1- exp(-TD./tauDD1 ) ) -rfi./alpha .* exp( -TD./tauDD1 ) ) ...
    ./UAD1 ...
    .*exp(-tpd./tau_U)...
    .*238/NA*1e9;

% U_D_nocc_extrema(2,2)=min(UBSED);
% U_D_nocc_extrema(2,4)=min(UBSED(Xei==3.2e7));
% U_D_nocc_extrema(2,6)=min(UBSED(Xei==3.2e8));
% U_D_nocc_extrema(2,8)=min(UBSED(n==n(3)));

for kn=1:Nn
    ntest=nset(kn);
    idx_all = find(n == ntest);
    %
    U_D_nocc_extrema(kn,15)=ntest;
    U_D_nocc_extrema(kn,16)=niso;    
    [U_D_nocc_extrema(kn,17),idx_sub]=min(UBSED(idx_all));
    index=idx_all(idx_sub);
    U_D_nocc_extrema(kn,18)=Xepd(index);U_D_nocc_extrema(kn,19)=Xei(index);U_D_nocc_extrema(kn,20)=rpd(index);
    U_D_nocc_extrema(kn,25)=min(UBSED(n==ntest & Xei==3.2e7));
    U_D_nocc_extrema(kn,27)=min(UBSED(n==ntest & Xei==3.2e8));
end

% max
rfpd=Xe_obs1.Xe132v130(2); % max fXe/130Xe

rfTD = ( ...
    UAD1.* exp(TD./tau_U)./ beta.* ( rfpd - rfa2.*( 1 - exp( -(tpd-TD)./tauDD2 ) ) ) ...
    +UAD2.* ( n.*rfa1.* ( 1 - exp( -TD./tauDD1 ) ) + rfi./alpha .* exp( -TD./tauDD1 )  )...
    )./ ( UAD2+ UAD1.* exp( TD./tau_U- (tpd-TD)./tauDD2 )./ beta );

UBSED=Xepd.* ...
    ( rfTD - n.*rfa1.*( 1- exp(-TD./tauDD1 ) ) -rfi./alpha .* exp( -TD./tauDD1 ) ) ...
    ./UAD1 ...
    .*exp(-tpd./tau_U)...
    .*238/NA*1e9;

% U_D_nocc_extrema(2,3)=max(UBSED);
% U_D_nocc_extrema(2,5)=max(UBSED(Xei==3.2e7));
% U_D_nocc_extrema(2,7)=max(UBSED(Xei==3.2e8));
% U_D_nocc_extrema(2,9)=max(UBSED(n==n(3)));

for kn=1:length(nset)
    ntest=nset(kn);
    idx_all = find(n == ntest);
    %   
    [U_D_nocc_extrema(kn,21),idx_sub]=max(UBSED(idx_all));
    index=idx_all(idx_sub);
    U_D_nocc_extrema(kn,22)=Xepd(index);U_D_nocc_extrema(kn,23)=Xei(index);U_D_nocc_extrema(kn,24)=rpd(index);
    U_D_nocc_extrema(kn,26)=max(UBSED(n==ntest & Xei==3.2e7));
    U_D_nocc_extrema(kn,28)=max(UBSED(n==ntest & Xei==3.2e8));
end

% 134
niso=134;
% U_D_nocc_extrema(3,1)=niso;
rfi=Xe134ren;
rfa=Xe134rpa;
rfa1=Xe134rpa;
rfa2=Xe134rpa;
YPu=YPu134;
YU=YUr134;

KPu1=YPu./(tau_Pu./tauDD1-1);
KU1=YU./(tau_U./tauDD1-1);
KPu2=YPu./(tau_Pu./tauDD2-1);
KU2=YU./(tau_U./tauDD2-1);

UAD1=KPu1 .* R0Pu_Ur .* ( exp(-TD./tau_Pu) - exp(-TD./tauDD1) ) ...
     +KU1 .* ( exp(-TD./tau_U) - exp(-TD./tauDD1) );
UAD2=KPu2 .* RTPu_Ur .* ( exp(-(tpd-TD)./tau_Pu) - exp(-(tpd-TD)./tauDD2) ) ...
     +KU2 .* ( exp(-(tpd-TD)./tau_U) - exp(-(tpd-TD)./tauDD2) );

% min
rfpd=Xe_obs1.Xe134v130(1); % min fXe/130Xe

rfTD = ( ...
    UAD1.* exp(TD./tau_U)./ beta.* ( rfpd - rfa2.*( 1 - exp( -(tpd-TD)./tauDD2 ) ) ) ...
    +UAD2.* ( n.*rfa1.* ( 1 - exp( -TD./tauDD1 ) ) + rfi./alpha .* exp( -TD./tauDD1 )  )...
    )./ ( UAD2+ UAD1.* exp( TD./tau_U- (tpd-TD)./tauDD2 )./ beta );

UBSED=Xepd.* ...
    ( rfTD - n.*rfa1.*( 1- exp(-TD./tauDD1 ) ) -rfi./alpha .* exp( -TD./tauDD1 ) ) ...
    ./UAD1 ...
    .*exp(-tpd./tau_U)...
    .*238/NA*1e9;

% U_D_nocc_extrema(3,2)=min(UBSED);
% U_D_nocc_extrema(3,4)=min(UBSED(Xei==3.2e7));
% U_D_nocc_extrema(3,6)=min(UBSED(Xei==3.2e8));
% U_D_nocc_extrema(3,8)=min(UBSED(n==n(3)));

for kn=1:Nn
    ntest=nset(kn);
    idx_all = find(n == ntest);
    %
    U_D_nocc_extrema(kn,29)=ntest;
    U_D_nocc_extrema(kn,30)=niso;    
    [U_D_nocc_extrema(kn,31),idx_sub]=min(UBSED(idx_all));
    index=idx_all(idx_sub);
    U_D_nocc_extrema(kn,32)=Xepd(index);U_D_nocc_extrema(kn,33)=Xei(index);U_D_nocc_extrema(kn,34)=rpd(index);
    U_D_nocc_extrema(kn,39)=min(UBSED(n==ntest & Xei==3.2e7));
    U_D_nocc_extrema(kn,41)=min(UBSED(n==ntest & Xei==3.2e8));
end

% max
rfpd=Xe_obs1.Xe134v130(2); % max fXe/130Xe

rfTD = ( ...
    UAD1.* exp(TD./tau_U)./ beta.* ( rfpd - rfa2.*( 1 - exp( -(tpd-TD)./tauDD2 ) ) ) ...
    +UAD2.* ( n.*rfa1.* ( 1 - exp( -TD./tauDD1 ) ) + rfi./alpha .* exp( -TD./tauDD1 )  )...
    )./ ( UAD2+ UAD1.* exp( TD./tau_U- (tpd-TD)./tauDD2 )./ beta );

UBSED=Xepd.* ...
    ( rfTD - n.*rfa1.*( 1- exp(-TD./tauDD1 ) ) -rfi./alpha .* exp( -TD./tauDD1 ) ) ...
    ./UAD1 ...
    .*exp(-tpd./tau_U)...
    .*238/NA*1e9;

% U_D_nocc_extrema(3,3)=max(UBSED);
% U_D_nocc_extrema(3,5)=max(UBSED(Xei==3.2e7));
% U_D_nocc_extrema(3,7)=max(UBSED(Xei==3.2e8));
% U_D_nocc_extrema(3,9)=max(UBSED(n==n(3)));

for kn=1:length(nset)
    ntest=nset(kn);
    idx_all = find(n == ntest);
    %   
    [U_D_nocc_extrema(kn,35),idx_sub]=max(UBSED(idx_all));
    index=idx_all(idx_sub);
    U_D_nocc_extrema(kn,36)=Xepd(index);U_D_nocc_extrema(kn,37)=Xei(index);U_D_nocc_extrema(kn,38)=rpd(index);
    U_D_nocc_extrema(kn,40)=max(UBSED(n==ntest & Xei==3.2e7));
    U_D_nocc_extrema(kn,42)=max(UBSED(n==ntest & Xei==3.2e8));
end


% 136
niso=136;
% U_D_nocc_extrema(4,1)=niso;
rfi=Xe136ren;
rfa=Xe136rpa;
rfa1=Xe136rpa;
rfa2=Xe136rpa;
YPu=YPu136;
YU=YUr136;


% %
% alpha=Xepd./Xei;
% rT=( (1-n).*ri + n.*(1-alpha).*ra ) ./ (1-n.*alpha);
% x1=(1./alpha-n)./(1-n);
% x2=(rT-ra)./(rpd-ra);
% 
% %
% TD=n.*tpd.*log(x1)./log(x1.^n.*x2);
% tauDD=tpd./log(x1.^n.*x2);
% tauDD1=n.*tauDD;
% tauDD2=tauDD;

KPu1=YPu./(tau_Pu./tauDD1-1);
KU1=YU./(tau_U./tauDD1-1);
KPu2=YPu./(tau_Pu./tauDD2-1);
KU2=YU./(tau_U./tauDD2-1);

% RTPu_Ur = R0Pu_Ur .* exp( -TD./tau_Pu + TD./tau_U );

UAD1=KPu1 .* R0Pu_Ur .* ( exp(-TD./tau_Pu) - exp(-TD./tauDD1) ) ...
     +KU1 .* ( exp(-TD./tau_U) - exp(-TD./tauDD1) );
UAD2=KPu2 .* RTPu_Ur .* ( exp(-(tpd-TD)./tau_Pu) - exp(-(tpd-TD)./tauDD2) ) ...
     +KU2 .* ( exp(-(tpd-TD)./tau_U) - exp(-(tpd-TD)./tauDD2) );

% min
rfpd=Xe_obs1.Xe136v130(1); % min fXe/130Xe

rfTD = ( ...
    UAD1.* exp(TD./tau_U)./ beta.* ( rfpd - rfa2.*( 1 - exp( -(tpd-TD)./tauDD2 ) ) ) ...
    +UAD2.* ( n.*rfa1.* ( 1 - exp( -TD./tauDD1 ) ) + rfi./alpha .* exp( -TD./tauDD1 )  )...
    )./ ( UAD2+ UAD1.* exp( TD./tau_U- (tpd-TD)./tauDD2 )./ beta );

UBSED=Xepd.* ...
    ( rfTD - n.*rfa1.*( 1- exp(-TD./tauDD1 ) ) -rfi./alpha .* exp( -TD./tauDD1 ) ) ...
    ./UAD1 ...
    .*exp(-tpd./tau_U)...
    .*238/NA*1e9;
% 
% U_D_nocc_extrema(4,2)=min(UBSED);
% U_D_nocc_extrema(4,4)=min(UBSED(Xei==3.2e7));
% U_D_nocc_extrema(4,6)=min(UBSED(Xei==3.2e8));
% U_D_nocc_extrema(4,8)=min(UBSED(n==n(3)));

for kn=1:Nn
    ntest=nset(kn);
    idx_all = find(n == ntest);
    %
    U_D_nocc_extrema(kn,43)=ntest;
    U_D_nocc_extrema(kn,44)=niso;    
    [U_D_nocc_extrema(kn,45),idx_sub]=min(UBSED(idx_all));
    index=idx_all(idx_sub);
    U_D_nocc_extrema(kn,46)=Xepd(index);U_D_nocc_extrema(kn,47)=Xei(index);U_D_nocc_extrema(kn,48)=rpd(index);
    U_D_nocc_extrema(kn,53)=min(UBSED(n==ntest & Xei==3.2e7));
    U_D_nocc_extrema(kn,55)=min(UBSED(n==ntest & Xei==3.2e8));
end


% max
rfpd=Xe_obs1.Xe136v130(2); % max fXe/130Xe

rfTD = ( ...
    UAD1.* exp(TD./tau_U)./ beta.* ( rfpd - rfa2.*( 1 - exp( -(tpd-TD)./tauDD2 ) ) ) ...
    +UAD2.* ( n.*rfa1.* ( 1 - exp( -TD./tauDD1 ) ) + rfi./alpha .* exp( -TD./tauDD1 )  )...
    )./ ( UAD2+ UAD1.* exp( TD./tau_U- (tpd-TD)./tauDD2 )./ beta );

UBSED=Xepd.* ...
    ( rfTD - n.*rfa1.*( 1- exp(-TD./tauDD1 ) ) -rfi./alpha .* exp( -TD./tauDD1 ) ) ...
    ./UAD1 ...
    .*exp(-tpd./tau_U)...
    .*238/NA*1e9;

% U_D_nocc_extrema(4,3)=max(UBSED);
% U_D_nocc_extrema(4,5)=max(UBSED(Xei==3.2e7));
% U_D_nocc_extrema(4,7)=max(UBSED(Xei==3.2e8));
% U_D_nocc_extrema(4,9)=max(UBSED(n==n(3)));

for kn=1:length(nset)
    ntest=nset(kn);
    idx_all = find(n == ntest);
    %   
    [U_D_nocc_extrema(kn,49),idx_sub]=max(UBSED(idx_all));
    index=idx_all(idx_sub);
    U_D_nocc_extrema(kn,50)=Xepd(index);U_D_nocc_extrema(kn,51)=Xei(index);U_D_nocc_extrema(kn,52)=rpd(index);
    U_D_nocc_extrema(kn,54)=max(UBSED(n==ntest & Xei==3.2e7));
    U_D_nocc_extrema(kn,56)=max(UBSED(n==ntest & Xei==3.2e8));
end





%% collect  extrema U data for piecewise D with given n
%
Xepd=input_test_step(:,1);
Xei=input_test_step(:,2);
rpd=input_test_step(:,3);
n=input_test_step(:,4);
m=input_test_step(:,5);
alpha=Xepd./Xei;


%
% Nn=length(0.01:0.0098:0.99);
Nn=length(0.01:0.049:0.99);
nset=m(1:Nn);
Nniso=4;

%

UBSED_m=zeros(Nn,15);


%
for kn=1:Nn

    ntest=nset(kn);
    idx_all = find(n == ntest);
    
    %
    UBSED_m(kn,1)=ntest;
    idx_iso=[131,132,134,136];

    idx_min=[3,17,31,45];
    [UBSED_m(kn,2),idx_sub]=max(U_D_nocc_extrema(kn,idx_min)); % max min
    index=idx_iso(idx_sub);index_min=idx_min(idx_sub);
    UBSED_m(kn,3)=index;
    UBSED_m(kn,4)=U_D_nocc_extrema(kn,index_min+1);UBSED_m(kn,5)=U_D_nocc_extrema(kn,index_min+2);UBSED_m(kn,6)=U_D_nocc_extrema(kn,index_min+3);

    idx_max=[7,21,35,49];
    [UBSED_m(kn,7),idx_sub]=min(U_D_nocc_extrema(kn,idx_max)); % min max
    index=idx_iso(idx_sub);index_max=idx_max(idx_sub);
    UBSED_m(kn,8)=index;
    UBSED_m(kn,9)=U_D_nocc_extrema(kn,index_max+1);UBSED_m(kn,10)=U_D_nocc_extrema(kn,index_max+2);UBSED_m(kn,11)=U_D_nocc_extrema(kn,index_max+3);

    UBSED_m(kn,12)=max(U_D_nocc_extrema(kn,[11,25,39,53]));
    UBSED_m(kn,13)=min(U_D_nocc_extrema(kn,[12,26,40,54]));
    UBSED_m(kn,14)=max(U_D_nocc_extrema(kn,[13,27,41,55]));
    UBSED_m(kn,15)=min(U_D_nocc_extrema(kn,[14,28,42,56]));

end


%%

figure;
set(gcf,'color','white','Units','centimeters','Position',[2 5 17.8 12]);

ax1=subplot(3,2,1);
contourf(TD_m(:,:,1),TD_m(:,:,2),TD_m(:,:,3),7);
xlabel('n');ylabel('m');
title('Min Transition time');
set(ax1,'FontUnits','points','FontSize',8,'FontWeight','bold',...
    'LineWidth',1.5,'LabelFontSizeMultiplier',1.0);
caxis([0 tpd]);
text(-0.126,1.01,'a','Units','normalized','FontWeight','bold','FontSize',8);

ax2=subplot(3,2,2);
contourf(TD_m(:,:,1),TD_m(:,:,2),TD_m(:,:,7),7);
xlabel('n');ylabel('m');
title('Max Transition time');
caxis([0 tpd]);c=colorbar;title(c,'Gyr');
set(ax2,'FontUnits','points','FontSize',8,'FontWeight','bold',...
    'LineWidth',1.5,'LabelFontSizeMultiplier',1.0);
text(-0.126,1.01,'b','Units','normalized','FontWeight','bold','FontSize',8);

ax3=subplot(3,2,3);
contourf(tauDD1_m(:,:,1),tauDD1_m(:,:,2),tauDD1_m(:,:,3),7);
title('Min Early Degassing Timescale');
xlabel('n');ylabel('m');
caxis([0 0.8]);
set(ax3,'FontUnits','points','FontSize',8,'FontWeight','bold',...
    'LineWidth',1.5,'LabelFontSizeMultiplier',1);
text(-0.126,1.01,'c','Units','normalized','FontWeight','bold','FontSize',8);

ax4=subplot(3,2,4);
contourf(tauDD1_m(:,:,1),tauDD1_m(:,:,2),tauDD1_m(:,:,7),7);
xlabel('n');ylabel('m');
title('Max Early Degassing Timescale');
caxis([0 0.8]);c=colorbar;title(c,'Gyr');
set(ax4,'FontUnits','points','FontSize',8,'FontWeight','bold',...
    'LineWidth',1.5,'LabelFontSizeMultiplier',1);
text(-0.126,1.01,'d','Units','normalized','FontWeight','bold','FontSize',8);

ax5=subplot(3,2,5);
contourf(tauDD2_m(:,:,1),tauDD2_m(:,:,2),tauDD2_m(:,:,3),7);
title('Min Recent Degassing Timescale');
xlabel('n');ylabel('m');
caxis([0 2]);
set(ax5,'FontUnits','points','FontSize',8,'FontWeight','bold',...
    'LineWidth',1.5,'LabelFontSizeMultiplier',1);
text(-0.126,1.01,'e','Units','normalized','FontWeight','bold','FontSize',8);

ax6=subplot(3,2,6);
contourf(tauDD2_m(:,:,1),tauDD2_m(:,:,2),tauDD2_m(:,:,7),7);
xlabel('n');ylabel('m');
title('Max Recent Degassing Timescale');
caxis([0 2]);c=colorbar;title(c,'Gyr');
set(ax6,'FontUnits','points','FontSize',8,'FontWeight','bold',...
    'LineWidth',1.5,'LabelFontSizeMultiplier',1);
text(-0.126,1.01,'f','Units','normalized','FontWeight','bold','FontSize',8);

% ── move subplot ──
all_ax = [ax1 ax2 ax3 ax4 ax5 ax6];
for i = 1:6
    pos = get(all_ax(i),'Position');
    pos(3) = pos(3) * 0.88;   % 
    set(all_ax(i),'Position',pos);
end
% —— put colorbar ——
cb1 = colorbar(ax2);
title(cb1,'Gyr');
set(cb1,'Position',[0.92, 0.70, 0.018, 0.22]);
set(cb1,'FontSize',8);
cb2 = colorbar(ax4);
title(cb2,'Gyr');
set(cb2,'Position',[0.92, 0.39, 0.018, 0.22]);
set(cb2,'FontSize',8);
cb3 = colorbar(ax6);
title(cb3,'Gyr');
set(cb3,'Position',[0.92, 0.08, 0.018, 0.22]);
set(cb3,'FontSize',8);

% exportgraphics(gcf, 'test_figS1s1_analytical_step_nm.pdf', ...
%     'ContentType','image', ...   % vector/image
%     'BackgroundColor','white', ... % 
%     'Resolution',600); 


%


figure;
set(gcf,'color','white','Units','centimeters','Position',[2 5 17.8 12]);


%
ax1=subplot(2,2,1);
plot(nset,UBSED_m(:,2),'LineWidth',2,'color','r');
hold on;
plot(nset,UBSED_m(:,7),'LineWidth',2,'color','r');
hold on;
plot(nset,UBSE0_m(1).*nset./nset,'LineWidth',2,'color','k');
hold on;
plot(nset,UBSE0_m(2).*nset./nset,'LineWidth',2,'color','k');
hold on;
xlabel('n');
ylabel('U_{BSE} (ppb)');
title('Xe(t_0)=3.2e7-3.2e8 atoms/g')
ylim([0 50]);
set(ax1,'FontUnits','points','FontSize',7,'FontWeight','bold',...
    'LineWidth',1.5,'LabelFontSizeMultiplier',1);
text(-0.126,1.01,'a','Units','normalized','FontWeight','bold','FontSize',8);

%
legend('piecewise D','', ...
    'Constant D','',...
     'Position',[0.6,0.7,0.2,0.09],'fontsize',10);

%
ax2=subplot(2,2,3);
plot(nset,UBSED_m(:,12),'LineWidth',2,'color','r');
hold on;
plot(nset,UBSED_m(:,13),'LineWidth',2,'color','r');
hold on;
plot(nset,UBSE0_m(3).*nset./nset,'LineWidth',2,'color','k');
hold on;
plot(nset,UBSE0_m(4).*nset./nset,'LineWidth',2,'color','k');
hold on;
xlabel('n');
ylabel('U_{BSE} (ppb)');
title('Xe(t_0)=3.2e7 atoms/g')
ylim([0 50]);
set(ax2,'FontUnits','points','FontSize',7,'FontWeight','bold',...
    'LineWidth',1.5,'LabelFontSizeMultiplier',1);
text(-0.126,1.01,'b','Units','normalized','FontWeight','bold','FontSize',8);


ax3=subplot(2,2,4);
plot(nset,UBSED_m(:,14),'LineWidth',2,'color','r');
hold on;
plot(nset,UBSED_m(:,15),'LineWidth',2,'color','r');
hold on;
plot(nset,UBSE0_m(5).*nset./nset,'LineWidth',2,'color','k');
hold on;
plot(nset,UBSE0_m(6).*nset./nset,'LineWidth',2,'color','k');
hold on;
% legend( ...
%     'piecewise D','',...
%     'piecewise R','','Location','best');
xlabel('n');
ylabel('U_{BSE} (ppb)');
title('Xe(t_0)=3.2e8 atoms/g')
ylim([0 50]);
set(ax3,'FontUnits','points','FontSize',7,'FontWeight','bold',...
    'LineWidth',1.5,'LabelFontSizeMultiplier',1);
text(-0.126,1.01,'c','Units','normalized','FontWeight','bold','FontSize',8);


% exportgraphics(gcf, 'test_figS1s1_analytical_step_U.pdf', ...
%     'ContentType','image', ...   % vector/image
%     'BackgroundColor','white', ... % 
%     'Resolution',600); 


