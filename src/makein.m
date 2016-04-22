function [tempfilename] = makein(inputfile)
% This function creates the input .mat file for the main function.
% The .mat file produced will contain all the needed variables to begin
% a new porous convection run. All units are SI, unless otherwise
% specified. This function does not use a previously computed T-P field,
% so calculates a hydrostatic initial pressure field from the temperature
% field.
%
% Timothy Crone (tjcrone@gmail.com)

% read input file
readinput(inputfile);

% time stepping 
%adaptivetime = 1; % set to unity for adaptive time stepping
if adaptivetime==1
    %nstep = 80000; % number of steps to take with adaptive time stepping
    t = zeros(1,nstep); % initialize t vector for adaptive time stepping
else
    stepsize = 1e5; % step size in seconds
    runtime = 3e7; % total run time in seconds (3e9 is about 100 years)
    t = 0:stepsize:runtime-stepsize; % create time vector built from stepsize and runtime
    nstep = length(t); % number of steps required in model run
end
nout = nstep/20; % number of steps to output (must be divisor of nstep)

% domain geometry
%nx = 20; % number of grid cells in x-direction (columns)
%nz = 50; % number of grid cells in z-direction (rows)
%d = 100; % grid cell size (uniform grid, meters)
x = linspace(d/2,(nx-1)*d,nx);
z = linspace(d/2,(nz-1)*d,nz);
[~,Z] = meshgrid(x,z);

% some constants
rhom = 2950; % rock or grain density (basalt)
cm = 1004; % rock heat capacity (basalt)
lamdam = 2; % rock thermal conductivity (basalt)
alpham = 2e-5; % rock thermal expansion coefficient (basalt)
phi = ones(nz,nx)*0.03; % porosity
g = 9.8; % gravitational constant

% initial temperature conditions
%Tcold = 0;
%Thot = 600;
%Tcut = 400;
%zreset = 1000;
%T = Z*(Tcut-Tcold)/(zreset)+Tcold;
%T(T>Tcut)=Tcut;
%T(Z>zreset)=Thot;
%T(T>Thot) = Thot; % make sure no values are above Thot
%T(T<Tcold) = Tcold; % make sure no values are below Tcold
%T = Z*(Thot-Tcold)/(nz*d)+Tcold; % vertical gradient
%T = T + 2*(rand(nz,nx)-0.5).*(Thot-Tcold)./100; % add some randomness to initial T
T = Z*0+Tcold;
T(Z>zreset)=Thot;

% do restart stuff here
% if restart, load the T. if not, set the initial temperature conditions





%initial permeability
%kon = 3e-15;
%koff = 1e-32; 
kx = ones(nz,nx)*kon;  % permeability in x-direction
kz = ones(nz,nx)*kon;  % permeability in z-direction
kx(Z>zreset) = koff;
kz(Z>zreset) = koff;

% define permeability function
%kfunc = 0; % set to unity if using a permeability function
%kcall = '[kx,kz] = thermalcracking(nx,nz,Z,kon,koff,g,T1);';

% temperature boundary conditions (0=Neumann 1=Dirichlet)
% first row/column is value, second is type
Tbt = [ones(1,nx)*Tcold; ones(1,nx)*1]; % Dirichlet cold
%Tbb = [ones(1,nx)*Thot; ones(1,nx)*1]; % Dirichlet hot
Tbb = [ones(1,nx)*0; ones(1,nx)*0]; % Neumann zero
Tbr = [ones(nz,1)*0 ones(nz,1)*0]; % Neumann zero
Tbl = [ones(nz,1)*0 ones(nz,1)*0]; % Neumann zero
%Tbl = [ones(nz,1)*Thot ones(nz,1)*1]; % Dirichlet hot

% top boundary conduction
% set this variable to unity to have conduction across the top boundary
topconduction = 1;

% load or globalize thermodynamic tables
global TT PP RHO CP BETA ALPHA
if isempty(TT)
   load('../hydrotables/hydrotab7.mat');
end

% calculate starting pressure field uxing calcinitp
Ptop = 20e6; % average seafloor pressure at top of domain
[P,Pbound,dPdzbound,rhobound] = calcinitp(nx,nz,T,Tbt,Tbb,Ptop,TT, ...
    PP,RHO,g,d);

% pressure boundary conditions (0=Neumann 1=Dirichlet)
Pbt = [ones(1,nx).*Ptop;ones(1,nx)*1]; % open
Pbb = [ones(1,nx).*0;ones(1,nx)*0]; % closed
Pbr = [ones(nz,1).*0 ones(nz,1)*0]; % closed
Pbl = [ones(nz,1).*0 ones(nz,1)*0]; % closed

% error-check nout
if mod(nstep,nout) ~= 0 || mod(nout,1) ~= 0
   error('Sorry, nout is not a divisor of nstep, or is not an integer.  Input file not written.');
end

% save the output to a temporary file
inoutdir = '~/research/crackingfronts/in_out/';
[status, tempfilename] = system(sprintf('mktemp %sinput.XXXXXX', inoutdir));
%tempfilename = [tempfilename(1:end-1), '.mat'];


% set infilename based on steady
%if steady==1
%  infilename = [fileprefix,'_steady_in.mat'];
%else
%  infilename = [fileprefix,'_cracky_in.mat'];
%end

% define full infile name and check to see if it already exists
%fullinfilename = ['../in_out/',infilename];
%nameexist = fopen(fullinfilename, 'r+');
%if nameexist ~= -1
%   sure = input(sprintf('\nAre you sure you want to overwrite the existing file, %s? (y/n) ', fullinfilename),'s');
%   if isempty(strfind(sure,'y')) || length(sure) ~= 1
%      fclose(nameexist);
%      error('Input file not written.');
%   end
%end

% save variables to an input .mat file
save(tempfilename(1:end-1),'adaptivetime','t','nstep','nout','nx','nz','d','cm','lamdam','phi', ...
   'rhom','kx','kz','g','T','P','Tbb','Tbl','Tbr','Tbt','Ptop','Pbt','Pbb','Pbl','Pbr', ...
   'alpham','rhobound','Pbound','topconduction','zreset','Thot', ...
   'kon','koff','Z','steady', '-v7.3');
%fprintf('\nInput file %s written.\n\n',fullinfilename);
