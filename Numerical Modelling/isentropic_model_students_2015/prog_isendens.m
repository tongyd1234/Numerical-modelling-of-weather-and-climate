function [snew] = prog_isendens(sold,snow,unow,dthetadt)

% global variables
% -----------------
global idbg nx nb nz dt dx dtdx dth idthdt

    if (idbg==1)
       fprintf('Prognostic step: Isentropic mass density ...\n');
    end %if
    
    % declare snew
    snew = zeros(nx+2*nb,nz);
    
    % *** Exercise 2.1 isentropic mass density ***
    % *** time step for isentropic mass density ***
    % *** edit here ***
    i=nb+(1:nx);
    k=1:nz;
    snew(i,k)=sold(i,k)-dtdx/2*(snow(i+1,k).*(unow(i+1,k)+unow(i+2,k)) - snow(i-1,k).*(unow(i,k)+unow(i-1,k)));
    if (idthdt==1 & nargin > 3)
    snew(i,k) = snew(i,k)- 0.5*dt./dth.*(dthetadt(i,k+1)+dthetadt(i,k))...
        *(snow(i,k+1)-snow(i,k-1));
    end
    %
    % *** Exercise 2.1 isentropic mass density ***
