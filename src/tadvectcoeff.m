function [B,D] = tadvectcoeff(nx,nz,d,qx,qz,rhof,cf, ...
    Tbb,Tbt,Tbr,Tbl,rhobt,rhobb,rhobr,rhobl,cfbt,cfbb, ...
    cfbr,cfbl);
%this function builds matrix and column vector coefficients that
%are needed to solve the implicit form of the temperature advection-
%diffusion equation.  B and D are coefficients such that in the
%absence of the diffusive term, dT/dt = BT+D, using an upwind scheme.
%more detail can be found starting on note page 402.

%build the thermodynamic coefficient beta and alpha
%beta = rhos.*cs.*(1-phi) + rhof.*cf.*phi; %beta as defined on p. 402
%tau = calctau(T,phi,rhof,cf,rhos,cs,drhodT2,dcfdT2);
alpha = rhof.*cf;

%build logical masks
eposmask = qx(:,2:end)>=0;
enegmask = qx(:,2:end)<0;
wposmask = qx(:,1:end-1)>=0;
wnegmask = qx(:,1:end-1)<0;
nposmask = qz(1:end-1,:)>=0;
nnegmask = qz(1:end-1,:)<0;
sposmask = qz(2:end,:)>=0;
snegmask = qz(2:end,:)<0;

%compute matrices with elements of B
eout = -alpha.*qx(:,2:end).*eposmask./d;
sout = -alpha.*qz(2:end,:).*sposmask./d;
wout = alpha.*qx(:,1:end-1).*wnegmask./d;
nout = alpha.*qz(1:end-1,:).*nnegmask./d;
allout = eout+sout+wout+nout;
ein = [-alpha(:,2:end).*qx(:,2:end-1) zeros(nz,1)].*enegmask./d;
sin = [-alpha(2:end,:).*qz(2:end-1,:);zeros(1,nx)].*snegmask./d;
win = [zeros(nz,1) alpha(:,1:end-1).*qx(:,2:end-1)].*wposmask./d;
nin = [zeros(1,nx);alpha(1:end-1,:).*qz(2:end-1,:)].*nposmask./d;

%reshape element matrices and fill B
B = spdiags(reshape(allout,nx*nz,1),0,nx*nz,nx*nz);
B(nz+1,1) = eps;
ein = [ein(:,end);reshape(ein(:,1:end-1),(nx-1)*nz,1)];
B = spdiags(ein,nz,B);
sin = reshape(sin,nx*nz,1);
B = spdiags([99999;sin(1:end-1)],1,B);
nin = reshape(nin,nx*nz,1);
B = spdiags([nin(2:end);99999],-1,B);
win = [reshape(win(:,2:end),(nx-1)*nz,1);win(:,1)];
B = spdiags(win,-nz,B);

%create D with boundary info
D = zeros(nz,nx);
%left(east)
D(:,1) = cfbl.*Tbl(:,1).*rhobl.*qx(:,1).*wposmask(:,1)./d;
%right(west)
D(:,end) = D(:,end) + cfbr.*Tbr(:,1).*rhobr.*-qx(:,end).* ...
    enegmask(:,end)./d;
%top(north)
D(1,:) = D(1,:) + cfbt.*Tbt(1,:).*rhobt.*qz(1,:).* ...
    nposmask(1,:)./d;
%bottom(south)
D(end,:) = D(end,:) + cfbb.*Tbb(1,:).*rhobb.*-qz(end,:).* ...
    snegmask(end,:)./d;

%and reshape D
D = reshape(D,nx*nz,1);
