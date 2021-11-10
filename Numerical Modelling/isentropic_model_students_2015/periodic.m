function [phi] = periodic(phi,nx,nb)
%  This subroutine makes the array phi(n1,n2) periodic. At the left
%  and right border the number of 'nb' points is overwritten. The
%  periodicity of this operation is 'nx'.

        i = 1:nb;
        phi(i,:) = phi(nx+i,:);
        phi(nb+nx+i,:) = phi(nb+i,:);

