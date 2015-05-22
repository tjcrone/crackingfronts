function [A,C] = Tdiffcoeff(nx,nz,d,lamdam,Tbr,Tbl,Tbb,Tbt,topconduction)
%this function builds matrix and column vector coefficients that
%are needed to solve the implicit form of the temperature advection-
%diffusion equation.  A and C are coefficients such that in the 
%absence of fluid flow, dT/dt = AT+C, using a second-order central
%differencing scheme.  more detail can be found starting on note
%page 402 or so. lamdam is assumed to be scalar.

%dummy variables for stiffness matrix construction
betax = ones(nz,nx+1).*lamdam./(d^2);
betaz = ones(nz+1,nx).*lamdam./(d^2);

%build element matrices to fill A
ealpha = [betax(:,2:end-1) zeros(nz,1)];
salpha = [betaz(2:end-1,:);zeros(1,nx)];
walpha = [zeros(nz,1) betax(:,2:end-1)];
nalpha = [zeros(1,nx);betaz(2:end-1,:)];
allalpha = ealpha+salpha+walpha+nalpha;

%reshape element matrices and fill A
A = spdiags(-reshape(allalpha,nx*nz,1),0,nx*nz,nx*nz);
ealpha = [ealpha(:,end);reshape(ealpha(:,1:end-1),(nx-1)*nz,1)];
A(1,nz+1) = eps;
A = spdiags(ealpha,nz,A);
salpha = reshape(salpha,nx*nz,1);
A(1,2) = eps;
A = spdiags([99999;salpha(1:end-1)],1,A);
nalpha = reshape(nalpha,nx*nz,1);
A(2,1) = eps;
A = spdiags([nalpha(2:end);99999],-1,A);
walpha = [reshape(walpha(:,2:end),(nx-1)*nz,1);walpha(:,1)];
A(nz+1,1) = eps;
A = spdiags(walpha,-nz,A);

% compute C (left and right side only for benchmarking)
C = zeros(nz,nx);
AC = C;
% left side
dloc = find(Tbl(:,2)==1); %dirichlet boundary condition locations
C(dloc,1) = (2*lamdam/(d^2))*Tbl(dloc,1);
AC(dloc,1) = -(2*lamdam/(d^2));
% right side
dloc = find(Tbr(:,2)==1);
C(dloc,end) = (2*lamdam/(d^2))*Tbr(dloc,1);
AC(dloc,end) = -(2*lamdam/(d^2));
% bottom
dloc = find(Tbb(2,:)==1);
C(end,dloc) = (2*lamdam/(d^2))*Tbb(1,dloc);
AC(end,dloc) = -(2*lamdam/(d^2));

% top conduction
if topconduction == 1
   dloc = find(Tbt(2,:)==1);
   C(1,dloc) = (2*lamdam/(d^2))*Tbt(1,dloc);
   AC(1,dloc) = -(2*lamdam/(d^2));
end

% reshape C
C = reshape(C,nx*nz,1);

%add AC to A
AC = reshape(AC,nx*nz,1);
Adiag0 = spdiags(A,0);
A(1,1)=eps;
A = spdiags(Adiag0+AC,0,A);
