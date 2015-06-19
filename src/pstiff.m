function [A,B,C] = pstiff(nx,nz,d,Se,rhof,rhobt,rhobb,rhobr, ...
    rhobl,qx,qz,kx,kz,mu,g,T,Pbt,Pbb,Pbr,Pbl,biot,dexdt,dezdt)
%this function builds matrix and column vector coefficients that
%are needed to solve the implicit form of the pressure equation.
%this formulation assumes no flow side boundaries, dirichlet top
%boundary and neumann bottom boundary.  this is fixed for now!

%load or globalize thermodynamic tables
global TT PP RHO CP BETA ALPHA
if isempty(TT)
    load('../hydrotables/hydrotab7.mat');
end

%Compute interface viscosities from mu
[mux,muz] = interfacemean(mu);

%Compute interface permeabilities
[kxx,dummy] = interfacemean(kx);
[dummy,kzz] = interfacemean(kz);

%compute betas
betameanx = kxx./mux;
betameanz = kzz./muz;

%calculate upwind density values
[rhofN,rhofS,rhofE,rhofW] = upwind(rhof,rhobt,rhobb,rhobr, ...
    rhobl,qx,qz);

%build element matrices to fill A (boundaries zeroed for A)
ealpha = [betameanx(:,2:end-1) zeros(nz,1)]./(d^2)./Se.*(rhofE./rhof);
salpha = [betameanz(2:end-1,:);zeros(1,nx)]./(d^2)./Se.*(rhofS./rhof);
walpha = [zeros(nz,1) betameanx(:,2:end-1)]./(d^2)./Se.*(rhofW./rhof);
nalpha = [zeros(1,nx);betameanz(2:end-1,:)]./(d^2)./Se.*(rhofN./rhof);
allalpha = ealpha+salpha+walpha+nalpha;

%reshape element matrices and fill A
A = spdiags(-reshape(allalpha,nx*nz,1),0,nx*nz,nx*nz);
A(nz+1,1) = eps;
ealpha = [ealpha(:,end);reshape(ealpha(:,1:end-1),(nx-1)*nz,1)];
A = spdiags(ealpha,nz,A);
walpha = [reshape(walpha(:,2:end),(nx-1)*nz,1);walpha(:,1)];
A = spdiags(walpha,-nz,A);
A(2,1) = eps;
salpha = reshape(salpha,nx*nz,1);
A = spdiags([99999;salpha(1:end-1)],1,A);
nalpha = reshape(nalpha,nx*nz,1);
A = spdiags([nalpha(2:end);99999],-1,A);


%compute hydrostatic vector (B)
%compute density at top using Pbt and T(1,:)
rhotop = interptim(PP,TT,RHO,Pbt(1,:)./100000,T(1,:));
rhofz = [rhotop ; (rhof(1:end-1,:)+rhof(2:end,:))/2 ; rhobb];
rhofz(end,:) = 0; %zero out for neumman bottom boundary
B = betameanz(1:end-1,:)./Se*g/d.*rhofz(1:end-1,:).*(rhofN./rhof) - ...
    betameanz(2:end,:)./Se*g/d.*rhofz(2:end,:).*(rhofS./rhof);
B = reshape(B,nx*nz,1);

%compute boundary vector (C) (***no flow sides***)
topC = zeros(nz,nx);
topC(1,:) = -2*kz(1,:)./mu(1,:)./Se(1,:)./d./d.*(rhofN(1,:)./rhof(1,:));
topA = spdiags(A,0)+reshape(topC,nx*nz,1);
A(1,1) = eps;
A = spdiags(topA,0,A);
topC(1,:) = topC(1,:).*-Pbt(1,:); %***dirichlet top assumed***

%bottom boundary must be all neumann
botC = zeros(nz,nx);
botC(end,:) = kz(end,:)./mu(end,:)./Se(end,:)./d.*Pbb(1,:).* ...
    (rhofS(end,:)./rhof(end,:));
%botA = spdiags(A,0)-reshape(botC,nx*nz,1); %no changes to A for neumann bottom boundary
%A(1,1) = eps;
%A = spdiags(botA,0,A);
%botC(end,:) = botC(end,:).*Pbb_d(1,:); %effectively dirichlet

C = reshape(botC+topC,nx*nz,1);

%compute strain-rate vector
%D = -reshape(biot./Se.*(dexdt+dezdt),nx*nz,1);

%at this point, dPdt = A*P + B + C;
