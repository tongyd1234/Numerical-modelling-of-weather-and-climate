function [lheat,qv,qc,qr,rainnc,rainncv,nc,nr] = seifert(u,dthetadt,t,pres,snew,qv,qc,...
                         qr,exn,zhtold,zhtnow,rainnc,rainncv,nc,nr)
                     
% ***********************************************
% Two-moment microphysical scheme (Seifert, 2001/2006)
% adapted from COSMO, Annette Miltenberger and Lukas Papritz (2012)
% ***********************************************

% global variables
% ------------------------------------------
global nz dx nxb dt r r_v cp dth iern idthdt

% define constants
%svp1 = 0.6112;
svp2 = 17.67;
svp3 = 29.65;
svpt0 = 273.15;
ep2 = r/r_v;
xlv = 2.5E06;

% store specific humidity
qv_ini = qv;
%qr_ini = qr;
%qc_ini = qc;
%nr_ini = nr;
%nc_ini = nc;

%annette
% constants
% ----------

rho0 = 1.225;
rho_w = 1000;          % Materialdichte von Fluessigwasser 
L_wd = 2.4*10^6;      % Verdampfungswaerme
K_T = 2.500.*10.^(-2);  % Waermeleitfaehigkeit
c_r = 1./2;

% characteristics of cloud droplet distribution (cloud_nue1mue1)
nu = 1;             % parameters describing assumed distribution 
x_max = 2.6*10^(-10); % maximal droplet mass 
x_min = 4.20*10^(-15); % minimal droplet mass

% characteristics of rain droplet distribution (rainULI)
rain_x_min = 2.6*10^(-10);     % minimale Teilchenmasse
rain_x_max = 3.*10^(-6);        % maximale Teilchenmasse
a_geo = 1.24.*10^(-1); % Koeff. Geometrie
b_geo = 0.333333;      % Koeff. Geometrie
a_ven = 0.780000;      % Koeff. Ventilation (PK)
b_ven = 0.308000;      % Koeff. Ventilation (PK)
rain_nu = 0;           % Breiteparameter der Verteilung

% parameters for autoconversion
k_c = 9.44*10^9;    % Long-Kernel
k_1 = 600;          % Parameter fuer Phi-Fkt. (autoconversion)
k_2 = 0.68;         % Parameter fuer Phi-Fkt. (autoconversion)
k_au = k_c./(20.*x_max).*(nu+2).*(nu+4)./(nu+1).^2;   % autoconversion constant

% parameters for accretion
k_3 = 5*10^(-4);      % Parameter fuer Phi-Fkt. (accretion)
k_r = 5.78;           % Parameter Kernel (accretion)

% parameter for rain selfcollection and break-up
k_sc = k_c.*(nu+2)./(nu+1);                         % selfcollection constant
k_rr = 4.33;
k_br = 1000;
D_br = 1.1.*10^-3;

% parameters for rain evaporation and sedimentation
rain_cmu0 = 6;
rain_cmu1 = 30;
rain_cmu2 = 10.^3;
rain_cmu3 = 1.1.*10.^(-3);
rain_cmu4 = 1;
rain_cmu5 = 2;
N_sc = 0.710;          % Schmidt-Zahl (PK)
n_f = 0.333;           % Exponent von N_sc im Vent-koeff.
m_f = 0.5;             % Exponent von N_re im Vent-koeff.
nu_l = 1.460.*10.^(-5); % Kinem. Visc. von Luft
aa = 9.65;
bb = 10.3;
cc = 600;
alf = 9.65;
bet = 10.3;
gamma = 600;
%annette

% transpose input fields
rainnc_tr = rainnc';
rainncv_tr = rainncv';
t_tr = t';

% reset rain rate to zero
rainncv_tr(:)=0.;

% compute density
% ------------------------
k=1:nz;
i=1:nxb;
pii = exn/cp;
rho(:,k) = snew(:,k).*dth./(zhtnow(:,k+1)-zhtnow(:,k));
t = 0.5*(pii(i,k+1).*repmat(t_tr(k+1),nxb,1)+pii(i,k).*repmat(t_tr(k),nxb,1));
p = 0.5.*(pres(i,k)+pres(i,k+1));%1E05* (0.5*(pii(i,k+1)+pii(i,k))^(1004/287));
gam = 2.5E06./(1004.*0.5.*(pii(i,k)+pii(i,k+1))); % L / (cp_d * exn/cp)

%annette
rrho_c = (rho0./rho);
rrho_04 = (rho0./rho).^0.5;
% annette

f5 = svp2*(svpt0 - svp3)*xlv/cp;

% lukas
% vertical wind 
i=2:(nxb-1);
k=1:nz;
dz_dx(1:nxb,1:nz) = 0. ;
dz_dx(i,k) = (zhtnow(i+1,k+1)+zhtnow(i+1,k) - zhtnow(i-1,k+1)-zhtnow(i-1,k)) ./ (4.*dx);
%w(i,k) = 0.5*(zhtnow(i,k) + zhtnow(i,k+1) - zhtold(i,k) - zhtold(i,k+1)) / dt + 0.5 .* (u(i+1,k)+u(i,k)) .* dz_dx(i,k); 
w(i,k) = 0.5 .* (u(i+1,k)+u(i,k)) .* dz_dx(i,k); 
if(idthdt==1)
    w(i,k) = w(i,k) + 0.5 .* (dthetadt(i,k)+dthetadt(i,k+1)) .* snew(i,k) ./ rho(i,k);
end %if
w(1,k) = 0;
w(nxb,k) = 0;
%fprintf('w : %g\n',max(max(w)))
% lukas

%fprintf('max w: %g',max(max(w)));

%annette
% nucleation
%-----------
% HUCM continental case (Texas CCN)
% N_ccn = 1260.*10^6;
% N_max = 3000.*10^6;
% N_min = 300.*10^6;
% S_max = 20;
% k_ccn = 0.308;

wcb_min = 0.1;%0.1;
scb_min = 0.0;

T_3 = 273.2; %triple point water
e_3 = 6.1078.*100; % saturation vapor pressure at triple point
A_w = 17.2693882; % constant f. saturation vapor pressure (water)
B_w = 35.86; % constant f. saturation vapor pressure (water)
e_ws_vec  = @(ta) e_3.*exp(A_w.*(ta-T_3)./(ta-B_w));

ssw = r_v.*rho.*qv.*t./e_ws_vec(t)-1.0;

qr = qr.*rho;
qv = qv.*rho;
qc = qc.*rho;
nr = nr.*rho;
nc = nc.*rho;

w_cb(1:nxb,1:nz) = 0;
for k = 2:nz
    ind = (w(:,k) > wcb_min & ssw(:,k) >= scb_min );%& ssw(:,k) > ssw(:,min(k-1,1)));
    w_cb(ind,k) = w(ind,k);
end

% parameter for exponential decrease of N_ccn with height:
z0_nccn = 4000.0;     % up to this height (m) constant unchanged value:
z1e_nccn = 2000.0;    % height interval at which N_ccn decreases by factor 1/e above z0_nccn:

% characteristics of different kinds of prototype CN: intermediate case
N_cn0 = 5000.*10^6;
etas  = 0.8;          % soluble fraction

%Look-up tables (for r2 = 0.03 mum, lsigs = 0.4)
wcb_ind =  [0   0.5   1.0   2.5   5.0];
ncn_ind =  [0     50     100     200     400     800     1600     3200     6400].*10^6;

ltab_nuc = [0      0       0       0       0       0        0        0        0  ;
    0     37.2    67.1   119.5   206.7   340.5    549.4    549.4    549.4;
    0     39.0    77.9   141.2   251.8   436.7    708.7   1117.7   1117.7;
    0     42.3    84.7   169.3   310.3   559.5    981.7   1611.6   2455.6;
    0     44.0    88.1   176.2   352.3   647.8   1173.0   2049.7   3315.6];
ltab_nuc = ltab_nuc.*10^6;

% % hard upper limit for number conc that eliminates also unrealistic high value
% % that would come from the dynamical core
nc(w_cb>0) = min(nc(w_cb>0),N_cn0);

% N_cn depends on height (to avoid strong in-cloud nucleation)
zml_k = 0.5.*(zhtnow(:,1:nz)+zhtnow(:,2:nz+1));
n_cn = N_cn0.*min(exp((z0_nccn-zml_k)./z1e_nccn), 1.0);  % exponential decrease with height

nccn = zeros(nxb,nz); 
lp = find(w_cb > 0); 
for jj = lp'
    locy = find(wcb_ind>w_cb(jj), 1 )-1;
    locx = find(ncn_ind>n_cn(jj), 1 )-1;

    if ( isempty(locy))
        locy = 1;
    end
    if ( isempty(locx))
        locx = 1;
    end
    %assb(w_cb(jj),'w_cb',n_cn(jj),'n_cn')

    if (locx < 9) && (locy<5)
        indx0 = (n_cn(jj)-ncn_ind(locx))./(ncn_ind(locx+1)-ncn_ind(locx));
        indx1 = 1-indx0;
        indy0 = (w_cb(jj)-wcb_ind(locy))./(wcb_ind(locy+1)-wcb_ind(locy));
        indy1 = 1-indy0;
        nccn(jj) = (indx0.*ltab_nuc(locy,locx)+indx1.*ltab_nuc(locy,locx+1)).*indy0+...
            (indx0.*ltab_nuc(locy+1,locx)+indx1.*ltab_nuc(locy+1,locx+1)).*indy1;
    elseif (locx<9) && (locy>=5)
        indx0 = (n_cn(jj)-ncn_ind(locx))./(ncn_ind(locx+1)-ncn_ind(locx));
        indx1 = 1-indx0;
        locy = min(locy,5);
        nccn(jj) = indx0.*ltab_nuc(locy,locx)+indx1.*ltab_nuc(locy,locx+1);
    elseif (locy<5) && (locx>=9)
        indy0 = (w_cb(jj)-wcb_ind(locy))./(wcb_ind(locy+1)-wcb_ind(locy));
        indy1 = 1-indy0;
        locx = min(locx,9);
        nccn(jj) = ltab_nuc(locy,locx).*indy0+ltab_nuc(locy+1,locy).*indy1;
    else
        locy = 5;
        locx = 9;
        nccn(jj) = ltab_nuc(locy,locx);
    end
end

% If n_cn is outside the range of the lookup table values, resulting
% NCCN are clipped to the margin values. For the case of these margin values
% beeing larger than n_cn (which happens sometimes, unfortunately), limit NCCN by n_cn:
nccn = min(nccn, n_cn);

nuc_n = etas.*nccn-nc;
nuc_n = max(nuc_n,0.0);
nuc_q = min(nuc_n.*x_min,qv);
nuc_q(nuc_q<0) = 0;
nuc_n = nuc_q./x_min;

nc = nc+nuc_n;
qc = qc + nuc_q;
qv = qv - nuc_q;
%nucn_max=max(nuc_n,nucn_max);

% autoconversion, accretion, selfcollection, break-up
%----------------------------------------------------
% autoconversion

if max(max(qc > 0))==1
    
    au(1:nxb,1:nz) = 0;
    sc(1:nxb,1:nz) = 0;
    
    ind = (qc>0);
    x_c = min(max(qc(ind)./nc(ind),x_min),x_max);
    au(ind) = k_au.*qc(ind).^2.*x_c.^2.*dt.*rrho_c(ind);
    
    if max(max(qc > 10^(-6)))>0
        ind1 = (qc>10^-6);
        tau = min(max(1-(qc(ind1)./(qc(ind1)+qr(ind1))),10^-25),0.9);
        phi = k_1.*tau.^k_2.*(1-tau.^k_2).^3;
        au(ind1) = au(ind1).*(1+phi./(1-tau).^2);
    end
    
    au(ind) = max(min(qc(ind),au(ind)),0);
    
    sc(ind) = k_sc.*qc(ind).^2.*dt.*rrho_c(ind);   % selfcollection cloud droplets
    
    nr_au = au(ind)./x_max;
    nc_au = min(nc(ind),sc(ind));
    
    qc = qc-au;
    qr = qr+au;
    nr(ind) = nr(ind)+nr_au;
    nc(ind) = nc(ind)-nc_au; 
end

% accretion
ac = zeros(nxb,nz);
if max(max(qc > 0 & qr >0))>0
    ind = (qc>0 & qr>0);
    tau = min(max(1-qc(ind)./(qc(ind)+qr(ind)),10^-25),1);
    phi = (tau./(tau+k_3)).^4;
    ac(ind) = k_r.*qc(ind).*qr(ind).*phi.*rrho_04(ind).*dt;
    ac = min(qc,ac);
    
    x_c = min(max(qc./nc,x_min),x_max);
    nc_ac = min(nc,ac./x_c);
    
    qr = qr+ac;
    qc = qc-ac;
    nc = nc-nc_ac;
    
end

% self-collection rain / breakup
if max(max(qr > 0))==1
    ind = (qr>0);
    x_r = min(max(qr(ind)./nr(ind),rain_x_min),rain_x_max);
    D_r = a_geo.*x_r.^b_geo;
    
    % selfcollection
    sc = k_rr.*nr(ind).*qr(ind).*rrho_04(ind).*dt;
    
    % breakup
    br = sc.*0;
    if max(max(D_r > 0.3*10^(-3)))>0
        ind1 = (D_r>0.3.*10^(-3));
        phi1 = k_br.*(D_r(ind1)-D_br)+1;
        br(ind1) = phi1.*sc(ind1);
    end
    
    nr_sc = min(nr(ind),sc-br);
    
    nr(ind) = nr(ind)-nr_sc;
end

nr(nr<0) = 0;
nc(nc<0) = 0;
qc(qc<0) = 0;
qr(qr<0) = 0;

if (iern == 1)

% evaporation of rain droplets
% -----------------------------
e_d = qv.*r_v.*t;
e_sw = e_ws_vec(t);
s_sw = e_d./e_sw - 1;

if max(max((s_sw < 0) + (qr > 0) + (qc < 10^(-9))))==3 
    % condition for the occurence of evaporation
    ind = (((s_sw < 0) + (qr > 0) + (qc < 10^(-9)))==3);
    
    eva_q = zeros(nxb,nz);
    eva_n = zeros(nxb,nz);
    
    d_vtp = 8.7602.*10.^(-5).*t(ind).^(1.81)./p(ind);
    g_d = 4.0.*pi./(L_wd.^2./(K_T.*r_v.*t(ind).^2)+r_v.*t(ind)./(d_vtp.*e_sw(ind)));
    
    x_r = qr(ind)./(nr(ind)+10^(-20));
    x_r = min(max(x_r,rain_x_min), rain_x_max);
    
    D_m = a_geo.*x_r.^b_geo;
    
    mue(D_m <= rain_cmu3) = rain_cmu0.*tanh((4.*rain_cmu2.*...
        (D_m(D_m <= rain_cmu3)-rain_cmu3)).^rain_cmu5)+rain_cmu4;
    mue(D_m > rain_cmu3) = rain_cmu1.*tanh((rain_cmu2.*...
       (D_m(D_m > rain_cmu3)-rain_cmu3)).^rain_cmu5)+rain_cmu4;
   
    mue = mue';
   
    lam = (pi./6.*rho_w.*(mue+3).*(mue+2).*(mue+1)./x_r).^(1/3);
    
    gfak = 1.357940435+mue.*(0.3033273220+mue.*(-0.1299313363.*10.^(-1) + ...
        mue.*(0.4002257774.*10.^(-3) -mue.*0.4856703981.*10.^(-5))));
    
    f_q = a_ven+b_ven.*N_sc.^n_f.*(aa./nu_l.*rrho_04(ind)).^m_f.*gfak./sqrt(lam).* ...
        (1-1/2.*(bb./aa).*(lam./(cc+lam)).^(mue+5/2)...
        -1/8.*(bb./aa).^2.*(lam./(2.*cc+lam)).^(mue+5/2) ...
        -1/16.*(bb./aa).^3.*(lam./(3.*cc+lam)).^(mue+5/2) ...
        -5/127.*(bb./aa).^4.*(lam./(4.*cc+lam)).^(mue+5/2));
    
    
    gamma_eva(gfak > 0) = gfak(gfak>0).*(1.1.*10.^(-3)./D_m).*exp(-0.2.*mue);
    gamma_eva(gfak < 0) = 1;
    gamma_eva = gamma_eva';
       
    eva_q(ind) = -g_d.*c_r.*nr(ind).*(mue+1)./lam.*f_q.*s_sw(ind).*dt;
    eva_n(ind) = gamma_eva.*eva_q(ind)./x_r;
        
    eva_q = max(eva_q,0);
    eva_n = max(eva_n,0);
    eva_q = min(eva_q,qr);
    eva_n = min(eva_n,nr);
    
    qv = qv + eva_q;
    qr = qr - eva_q;
    nr = nr - eva_n;
end

end %if
  
% conversion of mixing ratios to mass densities
% -------------------------------------------------------------------------
qr = qr./rho;
qv = qv./rho;
qc = qc./rho;
nr = nr./rho;
nc = nc./rho;
%annette

% saturation adjustment
% ----------------------
es = eswat1(t).*100;
qvs = ep2.*es./(p-es);

%saturation adjustment: condensation/evaporation 
produc = (qv-qvs)./(1+p./(p-es).*qvs.*f5./(t-svp3).^2);

produc = max(produc,-qc); % no evaporation if no cloud water
produc(nc<=0) = min(0,produc(nc<=0)); % no condensation if no cloud droplets
produc = min(qv,produc); % limit condensation to qv

qc = qc+produc;
qc(nc<=0) = 0;
nc(qc<=0) = 0;
qv = qv-produc;

% %Limit rain drop size
nr = max(nr,qr./rain_x_max);
nr = min(nr,qr./rain_x_min);
%nc = max(nc,qc./x_max);
nc = min(nc,5000.*10^6);
 
nr(nr<0) = 0;
nc(nc<0) = 0;
qc(qc<0) = 0;
qr(qr<0) = 0;

%annette
% sedimentation of rain droplets
% ------------------------------
dzmin = 10^10;
% density correction for fall velocities
rhocorr = (rho0./rho).^0.5;
adz = 1./(zhtnow(:,2:nz+1)-zhtnow(:,1:nz)); % reciprocal vertical grid
dzmin = min(1.0./adz,dzmin);

qr = qr.*rho;
nr = nr.*rho;

dt_sedi = min(dt,0.7.*dzmin./20.0);
nt_sedi = max(ceil(max(max(dt./dt_sedi))),1);
dt_sedi = dt./nt_sedi;

for i = 1:nt_sedi
    v_n_rain(1:nxb,1:nz) = 0;
    v_q_rain(1:nxb,1:nz) = 0;
    q_flux(1:nxb,1:nz) = 0;
    n_flux(1:nxb,1:nz) = 0;
     
    k = 1:nz;
    if (max(max((qr(:,k) > 10^(-20))))>0)
        ind = (qr(:,k)>10^(-20));
        
        x_r = zeros(nxb,nz);
        x_r(ind) = qr(ind)./nr(ind);
        x_r(ind) = min(max(x_r(ind),rain_x_min),rain_x_max);
        D_m = (6./(rho_w.*pi).*x_r).^(1./3);
        
        mue = zeros(nxb,nz);
        if max(max(qc >= 10^(-20) & qr>10^(-20)))>0
            mue(qc>=10^(-20) & qr>10^(-20)) = (rain_nu+1)./b_geo-1;
        end
        if max(max(D_m(ind) <= rain_cmu3 & qc(ind) <= 10^(-20)))>0
            ind1 = (((D_m <= rain_cmu3)+(qr>10^(-20))+(qc <= 10^(-20)))==3);
            mue(ind1) = rain_cmu0.*tanh((4.*rain_cmu2.*(D_m(ind1)-rain_cmu3)).^2)+rain_cmu4;
        end
        if max(max(D_m(ind) > rain_cmu3 & qc(ind) <= 10^(-20)))>0
            ind2 = (((D_m > rain_cmu3)+(qr>10^(-20))+(qc <= 10^(-20)))==3);
            mue(ind2) = rain_cmu1.*tanh((rain_cmu2.*(D_m(ind2)-rain_cmu3)).^2)+rain_cmu4;
        end
        
        D_r = (D_m.^3./((mue+3).*(mue+2).*(mue+1))).^(1/3);
        
        v_n = alf-bet./(1+gamma.*D_r).^(mue+1);
        v_q = alf-bet./(1+gamma.*D_r).^(mue+4);
        v_n = v_n.*rhocorr;
        v_q = v_q.*rhocorr;
        v_n = max(v_n,0.1);
        v_q = max(v_q,0.1);
        v_n = min(v_n,20);
        v_q = min(v_q,20);
        v_n_rain = -v_q;  % fall velocity
        v_q_rain = -v_q;  % fall velocity
    else
        v_n_rain(1:nxb,1:nz) = 0;
        v_q_rain(1:nxb,1:nz) = 0;
    end
        
    % lower boundary condition for fall velocity
    v_n_rain(:,1) = v_n_rain(:,2); 
    v_q_rain(:,1) = v_q_rain(:,2);
        
    for k = nz-1:-1:1    
      
    	v_nv = 0.5.*(v_n_rain(:,k+1)+v_n_rain(:,k));
        v_qv = 0.5.*(v_q_rain(:,k+1)+v_q_rain(:,k));
        % assuming v_nv, v_qv always_negative
        c_nv = -v_nv.*adz(:,k).*dt_sedi;
        c_qv = -v_qv.*adz(:,k).*dt_sedi;
    
        kk = k;
        s_nv(1:nxb,1) = 0.0;
        s_nv(c_nv<=1) = v_nv(c_nv<=1).*nr(c_nv<=1,k);
        if max(c_nv > 1)>0
        	cflag(1:nxb) = 0;
        	while (max(c_nv > 1)>0 && kk<nz-1) 
            	ind = (c_nv>1);
                cflag(ind) = 1;
                s_nv(ind) = s_nv(ind)+nr(ind,kk)./adz(ind,kk);
                c_nv(ind) = (c_nv(ind)-1).*adz(ind,kk+1)./adz(ind,kk);
                kk  = kk+1;
            end
            s_nv(cflag) = s_nv(cflag)+nr(cflag,kk)./adz(cflag,kk).*min(c_nv(cflag),1.0);
            s_nv(cflag) = -s_nv(cflag)./dt_sedi;
        end

        kk = k;
        s_qv(1:nxb,1) = 0.0;
        s_qv(c_qv<=1) = v_qv(c_qv<=1).*qr(c_qv<=1,k);
        if max(c_qv > 1)>0
            cflag(1:nxb) = 0;
            while (max(c_qv > 1) && kk<nz-1)   
                ind = (c_qv>1);
                cflag(ind) = 1;
                s_qv(ind) = s_qv(ind)+qr(ind,kk)./adz(ind,kk);
                c_qv(ind) = (c_qv(ind)-1).*adz(ind,kk+1)./adz(ind,kk);
                kk  = kk + 1;
            end
            s_qv(cflag) = s_qv(cflag) + qr(cflag,kk)./adz(cflag,kk).*min(c_qv(cflag),1.0);
            s_qv(cflag) = -s_qv(cflag)./dt_sedi;
        end

        % Flux-limiter to avoid negative values
        n_flux(:,k) = max(s_nv,n_flux(:,k+1)-nr(:,k)./(adz(:,k).*dt_sedi));
        q_flux(:,k) = max(s_qv,q_flux(:,k+1)-qr(:,k)./(adz(:,k).*dt_sedi));
    end  

    % uppper boundary condition
    n_flux(:,nz) = 0.0;
    q_flux(:,nz) = 0.0;

    k = 1:nz-1;
    nr(:,k) = nr(:,k)+(n_flux(:,k)-n_flux(:,k+1)).*adz(:,k).*dt_sedi;
    qr(:,k) = qr(:,k)+(q_flux(:,k)-q_flux(:,k+1)).*adz(:,k).*dt_sedi;
    
    rainncv = -q_flux(:,1)' .* 3600./ rho_w * 1000.;  % mm/h
    rainnc = rainnc - q_flux(:,1)' .* dt_sedi ./ rho_w * 1000.;   % mm        
end

% Sedimentation rates seem to be too slow for the sedimentation velocities
% of 9 m/s

qr = qr./rho;
nr = nr./rho;

qv(qv<0) = 0;
qc(qc<0) = 0;
nr(nr<0) = 0;
nc(nc<0) = 0;
% nc(qc < 10^(-20)) = min(nc(qc < 10^(-20)),qc(qc < 10^(-20))./x_min);
% nr(qr < 10^(-20)) = min(nr(qr < 10^(-20)),qr(qr < 10^(-20))./rain_x_min);
%annette

%finally update all variables (except temperature)
lheat = gam.*(qv_ini - qv);

% % for debugging
% fprintf('qc : %g\n',max(max(qc)))
% fprintf('qr : %g\n',max(max(qr)))
% fprintf('nc : %g\n',max(max(nc)))
% fprintf('nr : %g\n',max(max(nr)))
