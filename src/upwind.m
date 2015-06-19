function [N,S,E,W] = upwind(X,Xbt,Xbb,Xbr,Xbl,qx,qz)
%this function computes the upwind values for the property X
%for use in upwind finite difference schemes that need upwind
%values of certain properties. X is the property needed, while
%Xbt, Xbb, Xbr, and Xbl are that property's boundary values. qx
%and qz are the velocities on the grid cell faces. The
%x-direction is positive downward, and the z-direction is
%is positive to the right. To clarify, if the south boundary of
%grid cell "i" has a positive velocity, than S at "i" will be X
%at "i". If the south boundary of "i" is negative, than S at "i"
%will be X at "i+1". Boundary conditions must be dirichlet.

%build logical masks
Eposmask = qx(:,2:end)>=0;
Enegmask = qx(:,2:end)<0;
Wposmask = qx(:,1:end-1)>=0;
Wnegmask = qx(:,1:end-1)<0;
Nposmask = qz(1:end-1,:)>=0;
Nnegmask = qz(1:end-1,:)<0;
Sposmask = qz(2:end,:)>=0;
Snegmask = qz(2:end,:)<0;

%obtain "out" values for each face of each cell
Eout = Eposmask.*X;
Wout = Wnegmask.*X;
Nout = Nnegmask.*X;
Sout = Sposmask.*X;

%obtain "in" values for each face of each cell
Ein = Enegmask.*[X(:,2:end) Xbr(:,1)];
Win = Wposmask.*[Xbl(:,1) X(:,1:end-1)];
Nin = Nposmask.*[Xbt(1,:);X(1:end-1,:)];
Sin = Snegmask.*[X(2:end,:);Xbb(1,:)];

%combine "in" and "out" values
E = Eout+Ein;
W = Wout+Win;
N = Nout+Nin;
S = Sout+Sin;
