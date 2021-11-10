function [qvnew,qcnew,qrnew] = prog_moisture(unow,qvold,qcold,qrold,...
                              qvnow,qcnow,qrnow,qvnew,qcnew,qrnew,dthetadt)

% global variables
% ----------------------
global nz idbg nb nx dt dtdx idthdt dth  

    if (idbg==1)
       fprintf('Prognostic step: Moisture scalars ...\n');
    end %if

    % declare
    qvnew = zeros(nx+2*nb,nz);
    qcnew = zeros(nx+2*nb,nz);
    qrnew = zeros(nx+2*nb,nz);


    % *** Exercise 4.1 moisture advection ***
    % *** edit here *** 
    %
    i=nb+(1:nx);
    %k=1:nz;
    k=2:nz-1;
    qvnew(i,k)=qvold(i,k)-dtdx/2*(qvnow(i+1,k)-qvnow(i-1,k)).*(unow(i+1,k)+unow(i,k));
    qcnew(i,k)=qcold(i,k)-dtdx/2*(qcnow(i+1,k)-qcnow(i-1,k)).*(unow(i+1,k)+unow(i,k));
    qrnew(i,k)=qrold(i,k)-dtdx/2*(qrnow(i+1,k)-qrnow(i-1,k)).*(unow(i+1,k)+unow(i,k));
    %qvnew(i,k)=qvold(i,k)-dtdx/2*(qvnow(i+1,k)-qvnow(i-1,k)).*(unow(i+1,k)+unow(i,k))+dt(EP+G);
    %qcnew(i,k)=qcold(i,k)-dtdx/2*(qcnow(i+1,k)-qcnow(i-1,k)).*(unow(i+1,k)+unow(i,k))+dt(G-CC-AC);
    %qrnew(i,k)=qrold(i,k)-dtdx/2*(qrnow(i+1,k)-qrnow(i-1,k)).*(unow(i+1,k)+unow(i,k))+dt(CC+AC-EP);
    if (idthdt==1 & nargin > 3)
    qvnew(i,k) = qvnew(i,k)- 0.5*dt./dth.*(dthetadt(i,k+1)+dthetadt(i,k))...
        .*(qvnow(i,k+1)-qvnow(i,k-1));
    qcnew(i,k) = qcnew(i,k)- 0.5*dt./dth.*(dthetadt(i,k+1)+dthetadt(i,k))...
        .*(qcnow(i,k+1)-qcnow(i,k-1));
    qrnew(i,k) = qrnew(i,k)- 0.5*dt./dth.*(dthetadt(i,k+1)+dthetadt(i,k))...
        .*(qrnow(i,k+1)-qrnow(i,k-1));
    end
    
    % *** Exercise 4.1  ***

