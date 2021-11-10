function [phi] = relax(phi,nx,nb,phi1,phi2)

    % relaxation is done over nr grid points
    nr = 8;
    n=2*nb+nx;

    % initialize relaxation array
    % ---------------------------
    rel = [1 0.99 0.95 0.8 0.5 0.2 0.05 0.01]; % nr points
    
    % relaxation boundary conditions
    % ------------------------------
    if (ndims(phi) == 2)
                
                for i=1:nr
		        phi(i,:) = phi1*rel(i) + phi(i,:)*(1-rel(i));
		        phi(n+1-i,:) = phi2*rel(i) + phi(n+1-i,:)*(1-rel(i));
                end % for
    else
                
                for i=1:nr
                	phi(i) = phi1*rel(i) + phi(i)*(1-rel(i));
                	phi(n+1-i) = phi2*rel(i) + phi(n+1-i)*(1-rel(i));
                 end % for
    end %if





