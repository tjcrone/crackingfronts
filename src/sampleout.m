function [tout, Tout, crackedout, qzout] = sampleout()
% This function samples a set of output files in time.
%
% Timothy Crone

perm = '215e16';

files = {sprintf('../in_out/k%s_stead02_out', perm) ...
  sprintf('../in_out/k%s_uncrack01_out', perm) ...
  sprintf('../in_out/k%s_uncrack02_out', perm) ...
  sprintf('../in_out/k%s_uncrack03_out', perm) ...
  sprintf('../in_out/k%s_uncrack04_out', perm)};% ...
  %sprintf('../in_out/k%s_crack05_out', perm)};

% load first file
load(files{1});
Tout1 = Tout;
crackedout1 = crackedout;
qzout1 = qzout;
tout1 = tout;
[nz1, nx1, nt1] = size(Tout);

% load second file
load(files{2});
[nz2, nx2, nt2] = size(Tout);

% start building output matrices
Tall = ones(nz2,nx2,nt2+nt1-1)*Tout(end,1,1);
Tall(1:nz1,1:nx1,1:nt1) = Tout1;
Tall(:,:,nt1+1:end) = Tout(:,:,2:end);
crackedall = logical(Tall*0);
crackedall(1:nz1,1:nx1,1:nt1) = crackedout1;
crackedall(:,:,nt1+1:end) = crackedout(:,:,2:end);
qzall = zeros(nz2+1,nx2,nt2+nt1-1);
qzall(1:nz1+1,1:nx1,1:nt1) = qzout1;
qzall(:,:,nt1+1:end) = qzout(:,:,2:end);
tall = [tout1 tout(2:end)+tout1(end)];

% loop through remaining files
for i = 3:length(files)
  load(files{i});
  Tall = cat(3, Tall, Tout(:,:,2:end));
  crackedall = cat(3, crackedall, crackedout(:,:,2:end));
  qzall = cat(3, qzall, qzout(:,:,2:end));
  tall = [tall tout(2:end)+tall(end)];
end

% clear some variables to free up memory
clear Pout Tout Tout1 cfout crackedout crackedout1 qxout qzout qzout1 rhofout

% define interpolation times
tinterval = 100; % years
tmax = tall(end);
%if tmax/60/60/24/365/1000>325
%  tmax = 325*1000*365*24*60*60;
%end
tout = 0:tinterval*365*24*60*60:tmax;

% interpolate Tout
[X,Y,Z] = meshgrid(1:nx2,1:nz2,tall);
[X2,Y2,Z2] = meshgrid(1:nx2,1:nz2,tout);
Tout = interp3(X,Y,Z,Tall,X2,Y2,Z2);
clear Tall;

% interpolate crackedout
crackedout= interp3(X,Y,Z,crackedall,X2,Y2,Z2, 'nearest');
clear crackedall;

% interpolate qzout
[X,Y,Z] = meshgrid(1:nx2,1:nz2+1,tall);
[X2,Y2,Z2] = meshgrid(1:nx2,1:nz2+1,tout);
qzout= interp3(X,Y,Z,qzall,X2,Y2,Z2);
