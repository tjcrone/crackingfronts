function [intmeanx,intmeanz] = interfacemean(Y)
%This function computes the harmonic mean of any parameter
%(Y) on the control volume interfaces in the x and the z
%directions.  This function can be used to determine interface
%values for the viscosity, thermal conductivity and permeability
%in the porous convection codes.  Values on the boundary
%faces take on the value of the boundary nodes.  This function
%only works when the interfaces are midway between the nodes, but
%the grid does not have to be uniform.  See Patankar p. 44 for
%more on interface values and the harmonic mean.

warning off;
intmeanx = [Y(:,1) (2.*Y(:,1:end-1).*Y(:,2:end))./(Y(:,1:end-1) + ...
      Y(:,2:end)) Y(:,end)];
intmeanz = [Y(1,:); (2.*Y(1:end-1,:).*Y(2:end,:))./(Y(1:end-1,:) + ...
      Y(2:end,:));Y(end,:)];
warning on;

%where both node values are zero, this approach produces NANs for 
%the values on the interfaces.  we want zero on such an interface,
%so we set all NANs to zero
intmeanx(find(isnan(intmeanx)==1)) = 0;
intmeanz(find(isnan(intmeanz)==1)) = 0;