function [its_out] = ...
  makeoutput(unow,snow,zht,its_out,its,qvnow,qcnow,qrnow,ncnow,nrnow,dthetadt,tot_prec,prec)

% global variables
% ----------------
global dt nb nz nz1 nx nx1 nxb nxb1 r cp g idbg            ...
       imoist Z U S T QV QC QR NC NR TOT_PREC PREC LHEAT 

% compute
% ---------------
    if (idbg==1)
           fprintf('Prepare output ...\n');
    end %if
    
    k=1:nz;
    i=1:nx;
    
    uout(i,k) = 0.5*(unow(i+nb,k)+unow(i+nb+1,k));  % destagger
    sout(i,k) = snow(i+nb,k);

    if (imoist==1)
        qvout(i,k)=qvnow(i+nb,k);
        qcout(i,k)=qcnow(i+nb,k);
        qrout(i,k)=qrnow(i+nb,k);
        
        ncout(i,k)=ncnow(i+nb,k);
        nrout(i,k)=nrnow(i+nb,k);
    end %if
    
    if (imoist==1)
        lheatout(i,k)=0.5*(dthetadt(i+nb,k)+dthetadt(i+nb,k+1));    % destagger
    end %if
        
    if (imoist==1)
        rainout(i)    = prec(i+nb);
        totrainout(i) = tot_prec(i+nb);
    end %if
 
    if (idbg==1 || idbg==0)
        fprintf('Writing output ...\n');
     end %if
        
     % Prepare to write the 3-D field 
     % (2 spatial + 1 time) to logfile
     its_out = its_out+1 ;    
     Z(:,:,its_out) = zht(nb+1:nb+nx,:);
     U(:,:,its_out) = uout ;
     S(:,:,its_out) = sout ;
     if (imoist==1)
        QV(:,:,its_out) = qvout ;
        QC(:,:,its_out) = qcout ;
        QR(:,:,its_out) = qrout ;
        NC(:,:,its_out) = ncout ;
        NR(:,:,its_out) = nrout ;
        LHEAT(:,:,its_out) = lheatout;
        PREC(:,its_out) = rainout;
        TOT_PREC(:,its_out) = totrainout;
    end %if
    T(its_out) = its*dt ;



