function [initP,Pbound,dPdzbound,rhobound] = calcinitp(nx,nz,T,Tbt,Tbb,Ptop,TT,PP,RHO,g,d);
%this function calculates the initial pressure distribution within
%the model domain assuming a hydrostatic pressure gradient, and flow
%only in the z-direction.  at the moment, it only works for initial
%temperature distributions that are vertically uniform.  i expect
%that all initial temperature conditions will be vertically uniform,
%but if not, this function will have to be expanded. (see p. 388 for
%some derivations)  this function also computes the pressure gradient
%and the density of the inflowing fluid on the bottom boundary to be
%used if required

%initialize matrices
initP = [T*0+Ptop;zeros(1,nx)+Ptop];
coldP = initP;
hotP = initP;
Tcold = initP*0;
T = [T;Tbb(1,:)];
Pmid = initP;

%build and pressure A matrix
A = spdiags([-ones(nz+1,1),ones(nz+1,1)], -1:0, nz+1, nz+1);

%compute rhotop (cold)
rhotopcold = interptim(PP,TT,RHO,repmat(Ptop,1,nx)/1e5,T(1,:)*0);
%loop to find cold hydrostatic pressure (no flow)
for i = 1:100
   %compute rhof
   rhof = interptim(PP,TT,RHO,coldP./100000,Tcold);
   
   %build rhs   
   intrhof = (rhof(1:end-1,:)+rhof(2:end,:))/2;
   dPdz = intrhof*g*d;
   %dPdz = [(rhotopcold+rhof(1,:))/2*g*d/2+Ptop;dPdz];
   dPdz = [rhof(1,:)*g*d/2+Ptop;dPdz];
   dPdz(end,:) = intrhof(end,:)*g*d/2;
   
   %solve for P with backslash
   for j = 1:nx
       coldP(:,j) = A\dPdz(:,j);
   end
   %coldP = ([repmat(Ptop,1,nx);Pmid(1:end-1,:)]+Pmid)/2;
end

%compute rhotop (hot)
rhotophot = interptim(PP,TT,RHO,repmat(Ptop,1,nx)/1e5,T(1,:));
%loop to find hot hydrostatic (no flow)
for i = 1:100
    %compute rhof
   rhof = interptim(PP,TT,RHO,hotP./100000,T);
   
   %build rhs    
   intrhof = (rhof(1:end-1,:)+rhof(2:end,:))/2;
   dPdz = intrhof*g*d;
   %dPdz = [(rhotophot+rhof(1,:))/2*g*d/2+Ptop;dPdz];
   dPdz = [rhof(1,:)*g*d/2+Ptop;dPdz];
   dPdz(end,:) = intrhof(end,:)*g*d/2;
   
   %solve for P with backslash
   for j = 1:nx
       hotP(:,j) = A\dPdz(:,j);
   end
   %hotP = ([repmat(Ptop,1,nx);Pmid(1:end-1,:)]+Pmid)/2;
end

%calculate mean cold linearized hydrostatic gradient
alpha = repmat(mean(diff(coldP-hotP)/d),nz,1);

%loop to find hot hydrostatic + cold linear gradient
for i = 1:1000
    %compute rhof
   rhof = interptim(PP,TT,RHO,initP./100000,T);
   
   %build rhs
   intrhof = (rhof(1:end-1,:)+rhof(2:end,:))/2;
   dPdz = (intrhof*g+alpha)*d;
   %dPdz = [((rhotophot+rhof(1,:))/2*g+alpha(1,:))*d/2+Ptop;dPdz];
   dPdz = [(rhof(1,:)*g+alpha(1,:))*d/2+Ptop;dPdz];
   dPdz(end,:) = (intrhof(end,:)*g+alpha(1,:))*d/2;
   
   %solve for P with backslash
   for j = 1:nx
       initP(:,j) = A\dPdz(:,j);
   end
   %initP = ([repmat(Ptop,1,nx);Pmid(1:end-1,:)]+Pmid)/2;
end

%compute bottom boundary conditions (this assumes all boundaries
%are inflow boundaries, use the output from this section according
%to whether an inflow boundary is required
Pbound = initP(end,:);
rhobound = rhof(end,:);
dPdzbound = alpha(1,:);

%cut bottom value off of initP and rhof
initP = initP(1:end-1,:);

