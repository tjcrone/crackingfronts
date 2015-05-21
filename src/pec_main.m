function [] = pec_main(inputfile)
% This is the main poroelastic convection model function.  It loads the
% input variables specified by INPUTFILE, then commences a convection run
% saving the outputs to a file with the same prefix as the input file.
%
% Timothy Crone (tjcrone@gmail.com)

% load input file (no path, no extension)
load(['../in_out/',inputfile]);

% change time stepping for testing
%t = [0:1e6:1e10];
%nstep = length(t);
%nout = nstep;

% build outfile name
outfilename = inputfile(1:strfind(inputfile,'_')-1);
fulloutfilename = ['../in_out/',outfilename,'_fin'];

% add .mat for windows
if isempty(strfind(evalc('!uname'),'Linux'))
    fulloutfilename = [fulloutfilename,'.mat'];
end

% globalize thermodynamic tables
global TT PP RHO CP BETA
if isempty(TT)
    load('../hydrotables/hydrotab7.mat');
end

%add step signifier to some variables
T1 = T;
P1 = P;

%compute T-P dependent fluid properties (t=1)
mu1 = dynvisc(T1); %fluid viscosity
rhof1 = interptim(PP,TT,RHO,P1./100000,T1); %fluid density
cf1 = interptim(PP,TT,CP,P1./100000,T1); %fluid heat capacity
Se1 = T*0+1;

%compute boundary fluid properties
rhobt = interptim(PP,TT,RHO,Tbt(1,:)*0+Ptop./100000,Tbt(1,:));
rhobb = rhobound;
rhobr = interptim(PP,TT,RHO,P1(:,end)./100000,Tbr(:,1));
rhobl = interptim(PP,TT,RHO,P1(:,1)  ./100000,Tbl(:,1));
cfbt = interptim(PP,TT,CP,Tbt(1,:)*0+Ptop./100000,Tbt(1,:));
cfbb = interptim(PP,TT,CP,Pbound./100000,Tbb(1,:));
cfbr = interptim(PP,TT,CP,P1(:,end)./100000,Tbr(:,1));
cfbl = interptim(PP,TT,CP,P1(:,1)  ./100000,Tbl(:,1));

%compute darcy velocities (t=1)
[qx1,qz1] = pec_darcy(nx,nz,P1,rhof1,rhobb,kx,kz,mu1,g,d,Pbt,Pbb,Pbr,Pbl,T1);

%***Initialize Output Data Storage Locations***
rhofout = zeros(nz,nx,nout);
Tout = zeros(nz,nx,nout);
Pout = zeros(nz,nx,nout);
qxout = zeros(nz,nx+1,nout);
qzout = zeros(nz+1,nx,nout);
tout = zeros(1,nout);

%***Store t=1 Outputs
rhofout(:,:,1) = rhof1;
Tout(:,:,1) = T1;
Pout(:,:,1) = P1;
qxout(:,:,1) = qx1;
qzout(:,:,1) = qz1;

%create tentative values at t=2
rhof2 = rhof1;
cf2 = cf1;
qx2 = qx1;
qz2 = qz1;
Se2 = Se1;
mu2 = mu1;
P2 = P1;
T2 = T1;

% start timer
tic;

%time loop
for i = 1:nstep-1
    
    % calculate dt
    dt = t(i+1)-t(i);
   
    % compute beta for temperature equation
    beta2 = reshape(rhom.*cm.*(1-phi) + rhof2.*cf2.*phi,nx*nz,1);
        
    % create del-squared stiffness matrix for diffusion term in heat eq.
    [AimpT,CimpT] = Tdiffcoeff(nx,nz,d,lamdam,Tbr,Tbl,Tbb,Tbt);

    % compute T2 using implicit technique
    [BimpT,DimpT] = Tadvectcoeff2(nx,nz,d,qx2,qz2,rhof2,cf2, ...
        Tbb,Tbt,Tbr,Tbl,rhobt,rhobb,rhobr,rhobl,cfbt,cfbb, ...
        cfbr,cfbl);
        
    % single step implicit left hand side stiffness:
    Tstiff = spdiags(beta2,0,nx*nz,nx*nz) - dt*(AimpT + BimpT);
        
    % single step implicit right hand side
    RHS = reshape(T1,nx*nz,1).*beta2 + dt*(CimpT + DimpT); %see p. 464
        
    % and solve:
    T2 = Tstiff\RHS;
    T2 = reshape(T2,nz,nx);
    T2(T2<0) = 0; %a kluge to prevent negative temperatures

    % compute P2 using implicit technique
    [AimpP,BimpP,CimpP] = Pstiff(nx,nz,d,Se2,rhof2, ...
        rhobt,rhobb,rhobr,rhobl,qx2,qz2,kx,kz,mu2,g,T2,Pbt, ...
        Pbb,Pbr,Pbl);

    % single step implicit left hand side stiffness:
    Pstiffness = AimpP;

    % single step implicit right hand side
    PRHS = -BimpP-CimpP;

    % and solve:
    P2 = Pstiffness\PRHS;
    P2 = reshape(P2,nz,nx);
        
    % compute T-P dependent fluid properties (t=2)
    mu2 = dynvisc(T2);
    rhof2 = interptim(PP,TT,RHO,P2./100000,T2); %fluid density
    cf2 = interptim(PP,TT,CP,P2./100000,T2); %fluid heat capacity
        
    % compute darcy velocities (t=2)
    [qx2,qz2] = pec_darcy(nx,nz,P2,rhof2,rhobb,kx,kz,mu2,g,d,Pbt,Pbb,Pbr,Pbl,T2);

    %  





        
    % shift variables
    P1 = P2;
    T1 = T2;   
    
    % Store outputs
    if mod(i,nstep/nout) == 0;
        rhofout(:,:,i/(nstep/nout)+1) = rhof2;
        Tout(:,:,i/(nstep/nout)+1) = T2;
        Pout(:,:,i/(nstep/nout)+1) = P2;
        qxout(:,:,i/(nstep/nout)+1) = qx2;
        qzout(:,:,i/(nstep/nout)+1) = qz2;
        tout(i/(nstep/nout)+1) = t(i+1);
    end
    
    %Update progress bar
    progressbar(i,nstep-1,mfilename, 'working ...');
end

%***Print timing info to screen***
etime = toc;
fprintf('\n\nTotal time              %f seconds\n',etime);
fprintf('Number of model steps         %i steps\n',i);
fprintf('Time per model step     %f seconds\n\n',etime/nstep);

%save outputs to file
save(fulloutfilename,'rhofout','Tout','Pout','qxout','qzout','tout');
fprintf('\nOutput file %s written.\n\n',[outfilename,'_fin.mat']);
