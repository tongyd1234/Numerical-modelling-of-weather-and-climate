function [ncnew,nrnew] = prog_numdens(unow,ncold,nrold,...
                              ncnow,nrnow,ncnew,nrnew,dthetadt)

% global variables
% ----------------------
global nz idbg nb nx dt dx idthdt dth  

% local variables
% ----------------------
dtdx = dt/dx;

    if (idbg==1)
       fprintf('Prognostic step: Number densities ...');
    end %if

    % declare
    ncnew = zeros(nx+2*nb,nz);
    nrnew = zeros(nx+2*nb,nz);

    % *** Exercise 5.1/5.2 Two Moment Scheme ***
    % *** edit here *** 
    %
   i=nb+(1:nx);
   k=2:nz-1;
   ncnew(i,k)=ncold(i,k)-dtdx/2*(ncnow(i+1,k)-ncnow(i-1,k)).*(unow(i+1,k)+unow(i,k));
   nrnew(i,k)=nrold(i,k)-dtdx/2*(nrnow(i+1,k)-nrnow(i-1,k)).*(unow(i+1,k)+unow(i,k));
    
    if (idthdt==1 & nargin > 3)
    ncnew(i,k) = ncnew(i,k)- 0.5*dt./dth.*(dthetadt(i,k+1)+dthetadt(i,k))...
        .*(ncnow(i,k+1)-ncnow(i,k-1));
    nrnew(i,k) = nrnew(i,k)- 0.5*dt./dth.*(dthetadt(i,k+1)+dthetadt(i,k))...
        .*(nrnow(i,k+1)-nrnow(i,k-1));
    end
   
    %
    % *** Exercise 5.1/5.2  ***

    
