function [exn,mtg] = diag_montgomery(prs,mtg,th0,topo,topofact)

% global variables
% --------------------------------
global idbg pref cp nz nxb g dth rdcp 

    if (idbg==1)
       fprintf('Diagnostic step: Exner function and Montgomery potential ...\n');
    end %if
  
    % *** Exercise 2.2 (also 1.5) Diagnostic computation of Montgomery ***
    % *** Calculate Exner function and Montgomery potential ***
    % *** Edit here ***
    i=1:nxb;
    % add computation of Exner function:   
    
    exn=cp.*(prs./pref).^(rdcp);
    
    % add lower boundary condition at height mtg(i,1)
    mtg(i,1)=g*topo(i)*topofact+exn(i,1)*(dth/2+th0(1));

    % add integration loop upwards:
    for k=2:nz;
    mtg(i,k)=mtg(i,k-1)+dth*exn(i,k);
    end;
    % *** Exercise 2.2 Diagnostic computation  ***
