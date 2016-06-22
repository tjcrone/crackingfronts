function [tmpfilename] = makein(inputfile)
% This function creates the input .mat file for the main function from the simpler
% text input file. The .mat file produced will contain all the needed variables to
% begin a new porous convection run. All units are SI, unless otherwise specified.
%
% Timothy Crone (tjcrone@gmail.com)

% read input file
readinput(inputfile);

% time stepping 
if adaptivetime==1
    t = zeros(1,nstep); % initialize t vector for adaptive time stepping
else
    stepsize = 1e5; % step size in seconds
    runtime = 3e7; % total run time in seconds (3e9 is about 100 years)
    t = 0:stepsize:runtime-stepsize; % create time vector built from stepsize and runtime
    nstep = length(t); % number of steps required in model run
end
nout = nstep/outputinterval; % number of steps to output (must be divisor of nstep)

% check nout
if mod(nstep,nout) ~= 0 || mod(nout,1) ~= 0
   error('Output step number is not an integer. Check nstep and outputinterval.');
end

% domain geometry
x = linspace(d/2,(nx-1)*d,nx);
z = linspace(d/2,(nz-1)*d,nz);
[~,Z] = meshgrid(x,z);

% some constants
rhom = 2950; % rock or grain density (basalt)
cm = 1004; % rock heat capacity (basalt)
lamdam = 2; % rock thermal conductivity (basalt)
alpham = 2e-5; % rock thermal expansion coefficient (basalt)
phi = ones(nz,nx)*phi; % porosity
g = 9.8; % gravitational constant

% cracked boolean
cracked = logical(Z*0+1);
cracked(Z>frontdepth) = 0;

% initial temperature conditions
rng('default');
T = Z*0+Tcold+rand(nz,nx);
T(~cracked) = Thot;

% restart temperature and uncracked fields if required
if restart==1
  R = load(restartfile, 'Tout', 'crackedout', 'Pout');
  [m, n] = size(R.Tout(:,:,end));
  if n~=nx
    error('Restarted domains must have the same number of columns as the original.');
  end
  cracked(1:m,:) = R.crackedout(:,:,end);
  T(1:m,:) = R.Tout(:,:,end);
  % map Tres onto current geometry if necessary
  %if sum(size(T)==size(Tres))~=2
  %  if steady==0
  %    error('Geometry resizing only allowed when restarting into another steady state run.');
  %  end

  %if sum(size(T)==size(Tres))~=2 % changing geometry
  %  [m, n] = size(Tres);
  %  if n~=nx
  %    error('Restarted temperature field must have the same number of columns.');
  %  end
  %  T(1:m,:) = Tres;
  %else % same size as restart
  %  T = Tres;
  %  cracked(1:m,:) = R.crackedout(:,:,end);
  %end 

  %[m, n] = size(Tres);
  %if n~=nx
  %  error('Restarted temperature field must have the same number of columns.');
  %end
end

%initial permeability
kx = ones(nz,nx)*kon;  % permeability in x-direction
kz = ones(nz,nx)*kon;  % permeability in z-direction
kx(~cracked) = koff;
kz(~cracked) = koff;

% define permeability function
%kfunc = 0; % set to unity if using a permeability function
%kcall = '[kx,kz] = thermalcracking(nx,nz,Z,kon,koff,g,T1);';

% temperature boundary conditions (0=Neumann 1=Dirichlet)
% first row/column is value, second is type
Tbt = [ones(1,nx)*Tcold; ones(1,nx)*1]; % Dirichlet cold
Tbb = [ones(1,nx)*Tbottomvalue; ones(1,nx)*Tbottomtype];
Tbr = [ones(nz,1)*0 ones(nz,1)*0]; % Neumann zero
Tbl = [ones(nz,1)*0 ones(nz,1)*0]; % Neumann zero
%Tbl = [ones(nz,1)*Thot ones(nz,1)*1]; % Dirichlet hot

% load or globalize thermodynamic tables
global TT PP RHO CP
if isempty(TT)
   load('../hydrotables/hydrotab8.mat');
end

% calculate starting pressure field uxing calcinitp
%Ptop = 20e6; % average seafloor pressure at top of domain
[P,Pbound,dPdzbound,rhobound] = calcinitp(nx,nz,T,Tbt,Tbb,Ptop,TT, ...
    PP,RHO,g,d);
if restart==1
  P(1:m,:) = R.Pout(:,:,end);
  %if sum(size(P)==size(Pres))~=2
   % error('Geometry resizing not yet allowed');
  %else % same size as restart
  %  P = Pres;
  %end 
end

% pressure boundary conditions (0=Neumann 1=Dirichlet)
Pbt = [ones(1,nx).*Ptop;ones(1,nx)*1]; % open
Pbb = [ones(1,nx).*0;ones(1,nx)*0]; % closed
Pbr = [ones(nz,1).*0 ones(nz,1)*0]; % closed
Pbl = [ones(nz,1).*0 ones(nz,1)*0]; % closed

% save the output to a temporary file
inoutdir = '../in_out/';
[status, tmpfilename] = system(sprintf('mktemp %sinput.XXXXXX', inoutdir));

% save variables to an input .mat file
save(tmpfilename(1:end-1),'adaptivetime','t','nstep','nout','nx','nz','d','cm', ...
  'lamdam','phi','rhom','kx','kz','g','T','P','Tbb','Tbl','Tbr','Tbt','Ptop','Pbt', ...
  'Pbb','Pbl','Pbr','alpham','rhobound','Pbound','Ttopconduction','cracked','Thot', ...
  'kon','koff','Z','steady','maxpicard','picardthresh','tdamp','stepfraction', ...
  'Apress','Atemp','KIc','Pw','Tve','rhoR','-v7.3');
