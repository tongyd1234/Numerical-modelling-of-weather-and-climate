function write_logfile(nout)

global Z U S QV QC QR NR NC LHEAT TOT_PREC PREC T u00 bv00 th00 pref xl ...
       thl nx nz dx dth dt topomx topowd imoist     ...
       imicrophys idthdt imoist_diff run_name idbg

% Record output
% -----------------------------------
if (idbg==1 | idbg==0)
   fprintf('Writing to file %s \n',run_name);
   fprintf('Output contains %u output steps\n', nout )
end %if

logf.Z = Z ; % height in z coord.	
logf.U = U ; % velocity
logf.S = S ; % isentropic density
if (imoist==1)
logf.QV = QV ; % water vapor mix. ratio
logf.QC = QC ; % cloud mix. ratio
logf.QR = QR ; % rain mix. ratio
logf.LHEAT = LHEAT ; % latent heating
logf.PREC  = PREC ; % rain rate
logf.TOT_PREC  = TOT_PREC ; % accumulated rain
logf.NC = NC ; % number density of cloud droplets
logf.NR = NR ; % number density of rain droplets
end %if
logf.T = T ; % time

% variable parameters
v.u00 = u00;
v.bv00 = bv00;
v.th00 = th00;
v.pref = pref;
v.xl = xl;
v.thl = thl;
v.nx = nx;
v.nz = nz;
v.nt = nout;
v.dx = dx;
v.dth = dth;
v.dt = dt;
v.topomx = topomx;
v.topowd = topowd;

% namelist parameters
v.imoist       = imoist;
v.imoist_diff   = imoist_diff;
v.imicrophys = imicrophys;
v.idthdt   = idthdt;

logf.v = v;  
save(run_name,'logf');  
%save logtest.mat -ASCII U;

