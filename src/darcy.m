function [qx,qz] = darcy(nx,nz,P,rhof,rhobb,kx,kz,mu,g,d,Pbt,Pbb,Pbr,Pbl,T)
%This function computes the Darcy velocities in the x and z
%directions on the interfaces of control volumes in the model
%domain.  Model grid spacing must be uniform and specified by
%d.  For this version of the function, kx and kz can be non-uniform
%matrices.  g is expected to be a scalar.  P, rhof, and mu are 
%matrices, and Pb(x) are the pressure boundary conditions.

%load or globalize thermodynamic tables
global TT PP RHO CP BETA ALPHA
if isempty(TT)
   load('../hydrotables/hydrotab7.mat');
end

% compute interface viscosities from mu
[mux,muz] = interfacemean(mu);

% compute interface permeabilities
[kxx,dummy] = interfacemean(kx);
[dummy,kzz] = interfacemean(kz);

%Compute the internal interface darcy velocities
qx = -kxx(:,2:end-1)./mux(:,2:end-1) .* (P(:,2:end) - P(:,1:end-1)) ./ d;
qz = -kzz(2:end-1,:)./muz(2:end-1,:) .* ((P(2:end,:) - P(1:end-1,:))./d - g/2*(rhof(2:end,:)+rhof(1:end-1,:)));

%Pad qx and qz with dummy values to accept boundary velocities
qx = [zeros(nz,1) qx zeros(nz,1)];
qz = [zeros(1,nx);qz;zeros(1,nx)];

% Compute the boundary face darcy velocities using the boundary conditons
% Right Side:
nloc = find(Pbr(:,2)==0); %neumann boundary condition locations
dloc = find(Pbr(:,2)==1); %dirichlet boundary condition locations
floc = find(Pbr(:,2)==2); %flux boundary conditions
qx(nloc,end) = -kxx(nloc,end)./mux(nloc,end) .* (Pbr(nloc,1));
qx(dloc,end) = -kxx(dloc,end)./mux(dloc,end) .* (Pbr(dloc,1)-P(dloc,end))./(d/2);
qx(floc,end) = Pbr(floc,1);

% Left Side:
nloc = find(Pbl(:,2)==0); %neumann boundary condition locations
dloc = find(Pbl(:,2)==1); %dirichlet boundary condition locations
floc = find(Pbl(:,2)==2); %flux boundary conditions
qx(nloc,1) = -kxx(nloc,1)./mux(nloc,1) .* (Pbl(nloc,1));
qx(dloc,1) = -kxx(dloc,1)./mux(dloc,1) .* (P(dloc,1) - Pbl(dloc,1))./(d/2);
qx(floc,1) = Pbl(floc,1);

% Top:
nloc = find(Pbt(2,:)==0); %neumann boundary condition locations
dloc = find(Pbt(2,:)==1); %dirichlet boundary condition locations
floc = find(Pbt(2,:)==2); %flux boundary conditions
rhotop = zeros(1,nx);
rhotop(dloc) = interptim(PP,TT,RHO,Pbt(1,dloc)./100000,T(1,dloc));
qz(1,nloc) = -kzz(1,nloc)./muz(1,nloc) .* (Pbt(1,nloc)); %here, g is not subtracted because the non-hydrostatic pressure should have been input.
qz(1,dloc) = -kzz(1,dloc)./muz(1,dloc) .* ((P(1,dloc) - Pbt(1,dloc))./(d/2) - g*rhotop(dloc));
%qz(1,dloc) = -kzz(1,dloc)./muz(1,dloc) .* ((P(1,dloc) - Pbt(1,dloc))./(d/2) - g*rhof(1,dloc));
qz(1,floc) = Pbt(1,floc);

% Bottom:
nloc = find(Pbb(2,:)==0); %neumann boundary condition locations
dloc = find(Pbb(2,:)==1); %dirichlet boundary condition locations
floc = find(Pbb(2,:)==2); %flux boundary conditions
ndloc = find(Pbb(2,:)==3); %effectively dirichlet (ghost bin pressure)
qz(end,nloc) = -kzz(end,nloc)./muz(end,nloc) .* (Pbb(1,nloc)); %here, g is not subtracted because the non-hydrostatic pressure should have been input.
qz(end,dloc) = -kzz(end,dloc)./muz(end,dloc) .* ...
    ((Pbb(1,dloc) - P(end,dloc))./(d/2) - g*rhof(end,dloc));
qz(end,floc) = Pbb(1,floc);
qz(end,ndloc) = -kzz(end,ndloc)./muz(end,ndloc) .* ...
    ((Pbb(1,ndloc) - P(end,ndloc))./d - g*rhobb(ndloc));
