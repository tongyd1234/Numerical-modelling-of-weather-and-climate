function [eswat] = eswat1(T) 
%
% INPUT:  T     ... temperature [K]
%
% OUTPUT: Saturation vapor pressure over water [hPa]
%
%      Saturation vapor pressure over water, Goff-Gratch formulation
%      (based on exact integration of Clausius-Clapeyron equation)

%     define some constants:

      C1=7.90298;
      C2=5.02808;
      C3=1.3816E-7;
      C4=11.344;
      C5=8.1328E-3;
      C6=3.49149;


      RMIXV  = 373.16 ./ T;

      ES     = - C1*( RMIXV - 1) + C2*log10( RMIXV ) - ...
               C3*( 10.^( C4*(1 - 1./RMIXV)) - 1 ) +  ...
               C5*( 10.^( - C6*(RMIXV - 1)) - 1 );

      eswat  = 1013.246 * 10.^ES;
