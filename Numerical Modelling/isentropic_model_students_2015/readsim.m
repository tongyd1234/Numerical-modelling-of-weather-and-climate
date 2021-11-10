% -------------------------------------------------------------------
% Read the simulation data from the isentrop model
%
% Input:
%   simname		simulation name
%
% Output:
%   strucutre var with elements
%   .x[nx,nz]		time-independant x-coordinates [km]
%   .z[nx,nz,nt]	time-dependant z-coordinates [km]
%			(height of theta coordinate surfaces)
%   .u			horizontal velocity [m/s]
%   .w			vertical velocity [m/s]
%   .theta		theta [K]
%   .sigma		isentropic density
%   .rho		density
%   .topo		topography height
%   ...
%
% Example:
%   v = readsim('relax_nosponge');
%
% -------------------------------------------------------------------
function [v] = readsim(simname)
 
    load(simname); 

    % simulation parameters
    v = logf.v ;
    
    % display options
    v.xlim = [0, (v.nx-1)*v.dx/1000.]; 
    v.zlim = [0, 10];
    v.title = simname;
    v.tci = 2;		% theta contour interval (K)
    v.lhci = 1.0;	% latent heating (K/h) cont. interval
    v.vci = 2.;		% velocity contour interval (m/s)
    v.qvci = 0.5;	%vapor contour interval (g/kg)
    v.qcci = 0.1;	% cloud water contour interval (g/kg)
    v.qrci = 0.005;	% rain contour interval (g/kg)
    v.ncci = 10^8;	% cloud droplet number density contour interval (1/kg)
    v.nrci = 100;	% rain droplet number density contour interval (1/kg)
    v.tlim = [v.th00+v.dth/2, 400];        %K
    v.lhlim = [-5.0*v.lhci,5.0*v.lhci]; %K/h   contour limits
    v.vlim = [0,60];            %m/s
    v.qvlim = [v.qvci,10*v.qvci];      %g/kg
    v.qclim = [v.qcci,10*v.qcci];      %g/kg
    v.qrlim = [v.qrci,100*v.qrci];      %g/kg
    v.nclim = [v.ncci,10*v.ncci];	%1/kg
    v.nrlim = [v.nrci,10*v.nrci];	%1/kg
    v.totpreclim = [0.,0.5];    %mm
    v.preclim = [0.,1.];        %mm/h
    v.v0 = v.u00;
   
    % read output time step indices
    v.t = logf.T;

    % determine x-coordinate
    %v.x = repmat((1:v.nx)*v.dx, v.nz, 1)'; %debug
    v.x = repmat((0:v.nx-1)*v.dx/1000, v.nz, 1)'; %debug
    % z-coordinates (height of theta coordinate surfaces in km)
    v.zi = logf.Z/1000.;

    k = 1:v.nz;
    for i=1:v.nx
                %destagger 
		for t=1:v.nt
	    	v.z(i,k,t) = 0.5*(v.zi(i,k,t)+v.zi(i,k+1,t));
		end %for
    end %for
   
    % theta value between the coordinate surfaces
    v.theta = repmat(v.th00+v.dth/2:v.dth:(v.th00+v.thl),v.nx,1);

    % velocity field
    v.u = logf.U; %already destaggered

    % isentropic density
    v.sigma = logf.S;
    
    if (v.imoist==1) 
      v.qv    = logf.QV*1000; %convert to g/kg
      v.qc    = logf.QC*1000; %convert to g/kg
      v.qr    = logf.QR*1000; %convert to g/kg
      v.lheat = logf.LHEAT;   %(K/s)
      v.prec   =logf.PREC ;   % (mm/h)
      v.totprec=logf.TOT_PREC ;  % (mm)
      if (v.imicrophys==2)
        v.nc = logf.NC;	
        v.nr = logf.NR;
      end %if
    end %if

    % topography height
    v.topo = reshape(v.zi(:,1,:), [v.nx v.nt]); 

    % density
    k = 1:v.nz;
    for i=1:v.nx
	for t=1:v.nt
		v.rho(i,k,t) = v.sigma(i,k,t).*v.dth ... 
		   ./(v.zi(i,k+1,t)-v.zi(i,k,t));
	end %for
    end %for
    
    % axis ranges
    if (v.imoist == 1) 
        v.totpreclim = [0., max(ceil(max(max(v.totprec))*5.)/5., 0.2)];     %mm
    end %if