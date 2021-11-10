function [th0,exn0,prs0,z0,mtg0,s0,u0,sold,snow,uold,unow,mtg,mtgnew,...
          qv0,qc0,qr0,qvold,qvnow,qcold,qcnow,qrold,qrnow,...
          ncold,ncnow,nrold,nrnow]          =...
          makeprofile(sold,snow,uold,    ...
          unow,mtg,mtgnew,qvold,qvnow,qcold,qcnow,qrold,qrnow, ...
          ncold,ncnow,nrold,nrnow)
      
% Global variables
% ------------------------
global g cp cpdr pref th00 z00 exn00 u00 bv00    ...
idbg dth ishear k_shl k_sht u00_sh imoist imicrophys ...
nz nz1
k_c=12;
k_w=10;
rh_max=0.793;


if (idbg==1)
  fprintf('Create initial profile ...\n');
end %if


th0  = zeros(nz1,1);
exn0 = zeros(nz1,1);
bv0  = zeros(nz1,1);
prs0 = zeros(nz1,1);
z0   = zeros(nz1,1);
mtg0 = zeros(nz,1);
s0   = zeros(nz,1);
u0   = zeros(nz,1);
rh0  = zeros(nz,1);
qv0   = zeros(nz,1);
qc0   = zeros(nz,1);
qr0   = zeros(nz,1);
nc0   = zeros(nz,1);
nr0   = zeros(nz,1);


% Upstream profile for Brunt-Vaisala frequency (unstaggered)
% -------------------------------------------
k=1:nz1;
bv0(k)=bv00;

%  Upstream profile for theta (staggered)
% -------------------------------------------
k=1:nz1;
th0(k) = th00 + (k-1).*dth;

%  Upstream profile for Exner function and pressure (staggered)
% -------------------------------------------
exn0(1) = exn00;
for k=2:nz1
    exn0(k) = exn0(k-1) - 16 * g^2 * (th0(k)-th0(k-1)) / ((bv0(k-1) + bv0(k))^2 * (th0(k-1)+th0(k))^2);
end
k=1:nz1;
prs0(k) = pref .* (exn0(k)./cp).^cpdr;


% Upstream profile for geometric height (staggered)
% -------------------------------------------
z0(1) = z00;
for k=2:nz1
   z0(k) = z0(k-1) + 8*g*(th0(k)-th0(k-1))/((th0(k-1)+th0(k))*(bv0(k-1)+bv0(k))^2);
end


% Upstream profile for Montgomery potential (unstaggered)
% -------------------------------------------
mtg0(1) = g*z0(1) + th00*exn0(1) + dth*exn0(1)/2.;
k=2:nz;
mtg0(k) = mtg0(k-1) + dth.*exn0(k);

% Upstream profile for isentropic density (unstaggered)
% -------------------------------------------
k=1:nz;
s0(k) = -1./g .* (prs0(k+1)-prs0(k))./dth;

% Upstream profile for velocity (unstaggered)
% -------------------------------------------
k=1:nz;
u0(k) = u00;

% *** Exercise 3.3 Downslope windstorm ***
% *** use indices k_shl, k_sht, and wind speeds u00_sh, u00
%
if (ishear==1)
   if (idbg==1)
      fprintf('Using wind shear profile ...\n');
   end %if
% *** edit here ...


else
   if (idbg==1)
      fprintf('Using uniform wind profile ...\n');
   end %if
end %if
% *** Exercise 3.3 Downslope windstorm ***


% Upstream profile for moisture (unstaggered)
% -------------------------------------------
if (imoist ==1 && nargin>6)  % if number of function arguments > 6 -> moisture

% *** Exercise 4.1 Initial Moisture profile ***
% *** define new indices and create the profile ***
% *** for rh0; then use function rrmixv1 to compute qv0 ***
%
for k=1:nz
  if (k>(k_c-k_w) && k<(k_c+k_w))
    rh0(k)=rh_max*(cos(0.5*pi*abs(k-k_c)/k_w))^2;
  else 
      rh0(k)=0;
  end
    qv0(k)=rrmixv1(0.5.*(prs0(k)+prs0(k+1))./100,0.5.*(th0(k)./cp.*exn0(k)+th0(k+1)./cp.*exn0(k+1)),rh0(k),2);
end 
   
       

% *** edit here ...
% *** Exercise 4.1 Initial Moisture profile ***

end %if (nargin>6)

% Upstream profile for number densities (unstaggered)
% -------------------------------------------
if (imicrophys==2)
    nc0(k) = 0;
    nr0(k) = 0;
end %if



% Initial conditions for isentropic density (sigma), velocity u, and moisture qv
% ---------------------------------------------------------------------
for k=1:nz
    sold(:,k) = s0(k);
    snow(:,k) = s0(k);
    mtg(:,k) = mtg0(k);
    mtgnew(:,k) = mtg0(k);

    if (imoist==1 & nargin>6)
        qvold(:,k) = qv0(k);
        qvnow(:,k) = qv0(k);
        qcold(:,k) = qc0(k);
        qcnow(:,k) = qc0(k);
        qrold(:,k) = qr0(k);
        qrnow(:,k) = qr0(k);
    
        % droplet density for 2-moment scheme
        if (imicrophys==2)
             ncold(:,k) = nc0(k);
             ncnow(:,k) = nc0(k);
             nrold(:,k) = nr0(k);
             nrnow(:,k) = nr0(k);
        end % if
    end %if

    uold(:,k) = u0(k);
    unow(:,k) = u0(k);
end %for



