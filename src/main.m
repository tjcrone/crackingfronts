function [] = main(inputfile)
% This is the main poroelastic convection model function.  It loads the
% input variables specified by INPUTFILE, then commences a convection run
% saving the outputs to a file with the same prefix as the input file. If
% FINFILE is included, a restart is assumed and the final temperature from
% FINFILE is loaded as the starting temperature field.
%
% Timothy Crone (tjcrone@gmail.com)

% run makein
tempfilename = makein(inputfile);

% load the result of makein and delete the temporary mat file
load(tempfilename(1:end-1), '-mat');
system(sprintf('rm %s', tempfilename(1:end-1)));

% globalize thermodynamic tables
global TT PP RHO CP BETA
if isempty(TT)
    load('../hydrotables/hydrotab7.mat');
end

% add step signifier to some variables
T1 = T;
P1 = P;

% compute T-P dependent fluid properties (t=1)
mu1 = dynvisc(T1); %fluid viscosity
rhof1 = interptim(PP,TT,RHO,P1./100000,T1); %fluid density
cf1 = interptim(PP,TT,CP,P1./100000,T1); %fluid heat capacity
Se1 = T1*0+1;

% compute boundary fluid properties
rhobt = interptim(PP,TT,RHO,Tbt(1,:)*0+Ptop./100000,Tbt(1,:));
rhobb = rhobound;
rhobr = interptim(PP,TT,RHO,P1(:,end)./100000,Tbr(:,1));
rhobl = interptim(PP,TT,RHO,P1(:,1)  ./100000,Tbl(:,1));
cfbt = interptim(PP,TT,CP,Tbt(1,:)*0+Ptop./100000,Tbt(1,:));
cfbb = interptim(PP,TT,CP,Pbound./100000,Tbb(1,:));
cfbr = interptim(PP,TT,CP,P1(:,end)./100000,Tbr(:,1));
cfbl = interptim(PP,TT,CP,P1(:,1)  ./100000,Tbl(:,1));

% compute darcy velocities (t=1)
[qx1,qz1] = darcy(nx,nz,P1,rhof1,rhobb,kx,kz,mu1,g,d,Pbt,Pbb,Pbr,Pbl,T1);

% initialize output
rhofout = zeros(nz,nx,nout);
cfout = zeros(nz,nx,nout);
Tout = zeros(nz,nx,nout);
Pout = zeros(nz,nx,nout);
crackedout = zeros(nz,nx,nout);
qxout = zeros(nz,nx+1,nout);
qzout = zeros(nz+1,nx,nout);
tout = zeros(1,nout);

% store t=1 output
rhofout(:,:,1) = rhof1;
cfout(:,:,1) = cf1;
Tout(:,:,1) = T1;
Pout(:,:,1) = P1;
qxout(:,:,1) = qx1;
qzout(:,:,1) = qz1;
crackedout(:,:,1) = cracked;

% create tentative values at t=2
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

% time loop
for i = 1:nstep-1

    % if this is not a steady state run, crack
    if steady==0
      Apress=0.028;
      Atemp=4962.3; 
      KIc=1e6;
      Pw=200e5;
      Tve=900;
      rhoR=2900;
      KI=Atemp*(Tve-T1)-Apress*(Pw+rhoR*g*Z);
      cracked = ((KI>=KIc)+cracked)>0;
      kx = koff;
      kz = koff;
      kx(cracked) = kon;
      kz(cracked) = kon;
    end
    
    % set dt using adaptive or predefined time stepping
    if adaptivetime==1
        % adaptive time stepping based on CFL condition
        maxV = max(max(max(qx2)),max(max(qz2)));
        if i==1
            dt = 24*3600;
        elseif i<100
            dt = min(0.001*d/maxV, 0.5*d^2/1e-6);
        else
            dt = min(0.01*d/maxV, 0.5*d^2/1e-6);
        end
        dt = min([maxstepsize dt]);
        t(i+1)=t(i)+dt;
    else
        % use t vector for dt
        dt = t(i+1)-t(i);
    end

    % compute beta for temperature equation
    beta2 = reshape(rhom.*cm.*(1-phi) + rhof2.*cf2.*phi,nx*nz,1);
        
    % create del-squared stiffness matrix for diffusion term in heat eq.
    [AimpT,CimpT] = tdiffcoeff(nx,nz,d,lamdam,Tbr,Tbl,Tbb,Tbt,Ttopconduction);

    % compute T2 using implicit technique
    [BimpT,DimpT] = tadvectcoeff(nx,nz,d,qx2,qz2,rhof2,cf2, ...
        Tbb,Tbt,Tbr,Tbl,rhobt,rhobb,rhobr,rhobl,cfbt,cfbb, ...
        cfbr,cfbl);
        
    % single step implicit left hand side stiffness:
    Tstiff = spdiags(beta2,0,nx*nz,nx*nz) - dt*(AimpT + BimpT);
        
    % single step implicit right hand side
    RHS = reshape(T1,nx*nz,1).*beta2 + dt*(CimpT + DimpT); %see p. 464
        
    % and solve:
    T2 = Tstiff\RHS;
    T2 = reshape(T2,nz,nx);
    T2(T2<0) = 0; % kluge to prevent negative temperatures
    T2(T2>Thot) = Thot; % kluge to prevent overshoots

    % compute P2 using implicit technique
    [AimpP,BimpP,CimpP] = pstiff(nx,nz,d,Se2,rhof2, ...
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
    [qx2,qz2] = darcy(nx,nz,P2,rhof2,rhobb,kx,kz,mu2,g,d,Pbt,Pbb,Pbr,Pbl,T2);
        
    % shift variables
    P1 = P2;
    T1 = T2;   
    
    % store outputs
    if mod(i,nstep/nout) == 0;
        rhofout(:,:,i/(nstep/nout)+1) = rhof2;
        cfout(:,:,i/(nstep/nout)+1) = cf2;
        Tout(:,:,i/(nstep/nout)+1) = T2;
        Pout(:,:,i/(nstep/nout)+1) = P2;
        qxout(:,:,i/(nstep/nout)+1) = qx2;
        qzout(:,:,i/(nstep/nout)+1) = qz2;
        crackedout(:,:,i/(nstep/nout)+1) = cracked;
        tout(i/(nstep/nout)+1) = t(i+1);
    end
    
    % update progress bar
    progressbar(i,nstep-1,mfilename, 'working ...');
end

% print timing info
etime = toc;
fprintf('\n\nTotal wall time\t\t\t%.1f seconds\n',etime);
fprintf('Number of model steps\t\t%i steps\n',i+1);
fprintf('Wall time per step\t\t%.2f seconds\n',etime/nstep);
fprintf('Total model time\t\t%.1f years\n', tout(end)/60/60/24/365);
daysperstep = mean(diff(tout(end-round(length(tout)/4):end)))/60/60/24;
if daysperstep>365
    fprintf('Average model time per step\t%.2f years\n',daysperstep/365);
else
    fprintf('Average model time per step\t%.1f days\n',daysperstep);
end

% save outputs to file
underloc = strfind(inputfile, '_');
outfilename = [inputfile(1:underloc(end)-1), '_out.mat'];
save(outfilename,'rhofout', 'cfout', 'Tout','Pout','qxout','qzout','tout','crackedout','-v7.3');
fprintf('\nOutput file %s written.\n\n',outfilename);

