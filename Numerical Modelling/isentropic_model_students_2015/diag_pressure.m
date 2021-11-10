function [prs] = diag_pressure(prs0,prs,snow)

% global variables
% ---------------------------
global idbg nz nxb dth g

    if (idbg==1)
       fprintf('Diagnostic step: Pressure ...\n');
    end %if


    % *** Exercise 2.2 (also 1.4) Diagnostic computation of pressure ***
    % *** Diagnostic computation of pressure ***
    % *** (upper boundary condition and integration downwards) ***
    % *** edit here ***
    %
    prs(:,nz+1)=prs0(nz+1);
    i=1:nxb;
    for k=nz:-1:1;
    prs(i,k)=prs(i,k+1)+g*dth*snow(i,k);
    end;

    % *** Exercise 2.2 Diagnostic computation of pressure ***
