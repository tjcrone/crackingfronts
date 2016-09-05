function [tout, Tout, cfout, rhofout, crackedout, qzout] = sampleout()
% This function samples a set of output files in time.
%
% Timothy Crone

perm = '215e16';

files = {sprintf('../in_out/k%s_stead02_out', perm) ...
  sprintf('../in_out/k%s_crack01_out', perm) ...
  sprintf('../in_out/k%s_crack02_out', perm) ...
  sprintf('../in_out/k%s_crack03_out', perm) ...
  sprintf('../in_out/k%s_crack04_out', perm)};% ...
  %sprintf('../in_out/k%s_crack05_out', perm)};

% load first file
load(files{1});
Tout1 = Tout;
cfout1 = cfout;
rhofout1 = rhofout;
crackedout1 = crackedout;
qzout1 = qzout;
tout1 = tout;
t0 = tout(end);
tout1 = tout-t0;
[nz1, nx1, nt1] = size(Tout);

% load second file
load(files{2});
[nz2, nx2, nt2] = size(Tout);

% start building output matrices
Tall = ones(nz2,nx2,nt2+nt1-1)*Tout(end,1,1);
Tall(1:nz1,1:nx1,1:nt1) = Tout1;
Tall(:,:,nt1+1:end) = Tout(:,:,2:end);
cfall = ones(nz2,nx2,nt2+nt1-1)*cfout(end,1,1);
cfall(1:nz1,1:nx1,1:nt1) = cfout1;
cfall(:,:,nt1+1:end) = cfout(:,:,2:end);
rhofall = ones(nz2,nx2,nt2+nt1-1)*rhofout(end,1,1);
rhofall(1:nz1,1:nx1,1:nt1) = rhofout1;
rhofall(:,:,nt1+1:end) = rhofout(:,:,2:end);
crackedall = logical(Tall*0);
crackedall(1:nz1,1:nx1,1:nt1) = crackedout1;
crackedall(:,:,nt1+1:end) = crackedout(:,:,2:end);
qzall = zeros(nz2+1,nx2,nt2+nt1-1);
qzall(1:nz1+1,1:nx1,1:nt1) = qzout1;
qzall(:,:,nt1+1:end) = qzout(:,:,2:end);
tall = [tout1 tout(2:end)];

% loop through remaining files
for i = 3:length(files)
  load(files{i});
  Tall = cat(3, Tall, Tout(:,:,2:end));
  cfall = cat(3, cfall, cfout(:,:,2:end));
  rhofall = cat(3, rhofall, rhofout(:,:,2:end));
  crackedall = cat(3, crackedall, crackedout(:,:,2:end));
  qzall = cat(3, qzall, qzout(:,:,2:end));
  tall = [tall tout(2:end)+tall(end)];
end

% clear some variables to free up memory
%clear Pout Tout Tout1 cfout crackedout crackedout1 qxout qzout qzout1 rhofout

% define interpolation times
tinterval = 100; % years
tout1 = fliplr(-1*(0:tinterval*365*24*60*60:-min(tall)));
tout2 = 0:tinterval*365*24*60*60:max(tall);
tout = [tout1(1:end-1) tout2];
%tout = 0:tinterval*365*24*60*60:tmax;

% interpolate Tout
disp('Interpolating Tout.');
[X,Y,Z] = meshgrid(1:nx2,1:nz2,tall);
[X2,Y2,Z2] = meshgrid(1:nx2,1:nz2,tout);
Tout = interp3(X,Y,Z,Tall,X2,Y2,Z2);
%clear Tall;

% interpolate cfout
disp('Interpolating cfout.');
cfout= interp3(X,Y,Z,cfall,X2,Y2,Z2);

% interpolate rhofout
disp('Interpolating rhofout.');
rhofout= interp3(X,Y,Z,rhofall,X2,Y2,Z2);

% interpolate crackedout
disp('Interpolating crackedout.');
crackedout= interp3(X,Y,Z,crackedall,X2,Y2,Z2, 'nearest');
%clear crackedall;

% interpolate qzout
disp('Interpolating qzout.');
[X,Y,Z] = meshgrid(1:nx2,1:nz2+1,tall);
[X2,Y2,Z2] = meshgrid(1:nx2,1:nz2+1,tout);
qzout= interp3(X,Y,Z,qzall,X2,Y2,Z2);
