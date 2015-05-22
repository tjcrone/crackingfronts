function [] = makein()
% This function creates the input .mat file for the main function.
% The .mat file produced will contain all the needed variables to begin
% a new porous convection run. All units are SI, unless otherwise
% specified. This function does not use a previously computed T-P field,
% so computes a hydrostatic initial pressure field from the temperature
% field.
%
% Timothy Crone (tjcrone@gmail.com)

% infile name
infilename = 'testing04';

% time stepping variables
stepsize = 1e6; % step size in seconds
runtime = 1e9; % total run time in seconds (3e9 is about 100 years)
t = 0:stepsize:runtime-stepsize; % time vector built from stepsize and runtime
nstep = length(t); % number of steps required in model run
nout = nstep/10; % number of steps to output (must be divisor of nstep)

% domain geometry
nx = 25; % number of grid cells in x-direction (columns)
nz = 25; % number of grid cells in z-direction (rows)
d = 40; % grid cell size (uniform grid, meters)

% some constants
rhom = 2950; % rock or grain density (basalt)
cm = 1004; % rock heat capacity (basalt)
lamdam = 2; % rock thermal conductivity (basalt)
alpham = 2e-5; % rock thermal expansion coefficient (basalt)
phi = ones(nz,nx)*0.03; % porosity
kx = ones(nz,nx)*1e-12;  % permeability in x-direction
kz = ones(nz,nx)*1e-12;  % permeability in z-direction
g = 9.8; % gravitational constant

% initial temperature conditions
Tcold = 0;
Thot = 300;
x = linspace(d/2,(nx-1)*d,nx);
z = linspace(d/2,(nz-1)*d,nz);
[X,Z] = meshgrid(x,z);
T = Z*(Thot-Tcold)/(nz*d)+Tcold;
T = T + 2*(rand(nz,nx)-0.5).*(Thot-Tcold)./100; % add some randomness to initial T
T(T>Thot) = Thot; % make sure no values are above Thot
T(T<Tcold) = Tcold; % make sure no values are below Tcold

% temperature boundary conditions (0=Neumann 1=Dirichlet)
% first row/column is value, second is type
Tbt = [ones(1,nx)*Tcold; ones(1,nx)*1]; % Dirichlet cold
Tbb = [ones(1,nx)*Thot; ones(1,nx)*1]; % Dirichlet hot
%Tbb = [ones(1,nx)*0; ones(1,nx)*0]; % Neumann zero
Tbr = [ones(nz,1)*0 ones(nz,1)*0]; % Neumann zero
Tbl = [ones(nz,1)*0 ones(nz,1)*0]; % Neumann zero

% top boundary conduction
% set this variable to unity to have the conduction across the top boundary
topconduction = 0;

% load or globalize thermodynamic tables
global TT PP RHO CP BETA ALPHA
if isempty(TT)
   load('../hydrotables/hydrotab7.mat');
end

% calculate starting pressure field uxing calcinitp
Ptop = 20e6; % average seafloor pressure at top of domain
[initP,Pbound,dPdzbound,rhobound] = calcinitp(nx,nz,T,Tbt,Tbb,Ptop,TT, ...
    PP,RHO,g,d);
P = initP;
Pref = P;

% pressure boundary conditions (0=Neumann 1=Dirichlet)
Pbt = [ones(1,nx).*Ptop;ones(1,nx)*1]; % open
Pbb = [ones(1,nx).*0;ones(1,nx)*0]; % closed
Pbr = [ones(nz,1).*0 ones(nz,1)*0]; % closed
Pbl = [ones(nz,1).*0 ones(nz,1)*0]; % closed

% error-check nout
if mod(nstep,nout) ~= 0 | mod(nout,1) ~= 0
   error('Sorry, nout is not a divisor of nstep, or is not an integer.  Input file not written.');
end

% define full infile name and check to see if it already exists
fullinfilename = ['../in_out/',infilename,'_in'];
nameexist = fopen([fullinfilename,'.mat'],'r+');
if nameexist ~= -1
   sure = input(sprintf('\nAre you sure you want to overwrite the existing file, %s? (y/n) ', ...
      [infilename,'_in.mat']),'s');
   if isempty(findstr(sure,'y')) | length(sure) ~= 1
      fclose(nameexist);
      error('Input file not written.');
   end
end

% save variables to an input .mat file
save(fullinfilename,'stepsize','runtime','nstep','nout','nx','nz','d','cm','lamdam','phi', ...
   'rhom','kx','kz','g','T','P','Tbb','Tbl','Tbr','Tbt','Ptop','Pbt','Pbb','Pbl','Pbr', ...
   'alpham','Pref','rhobound','Pbound','t','topconduction');
disp(sprintf('\nInput file %s written.\n\n',[infilename,'_in.mat']));
