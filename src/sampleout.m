function [tI, TI, crackedI, qzI] = sampleout()
% This function samples a set of output files in time.
%
% Timothy Crone

files = {'../in_out/k215e15_stead02_out' ...
  '../in_out/k215e15_crack01_out' ...
  '../in_out/k215e15_crack02_out' ...
  '../in_out/k215e15_crack03_out' ...
  '../in_out/k215e15_crack04_out' ...
  '../in_out/k215e15_crack05_out'};

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
crackedall = Tall*0;
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

% interpolate
tinterval = 100; % years
tvec = 0:tinterval*365*24*60*60:tall(end);
[X,Y,Z] = meshgrid(1:nx2,1:nz2,tall);
[X2,Y2,tI] = meshgrid(1:nx2,1:nz2,tvec);
TI = interp3(X,Y,Z,Tall,X2,Y2,tI);
crackedI= interp3(X,Y,Z,crackedall,X2,Y2,tI);

[X,Y,Z] = meshgrid(1:nx2,1:nz2+1,tall);
[X2,Y2,tI] = meshgrid(1:nx2,1:nz2+1,tvec);
qzI= interp3(X,Y,Z,qzall,X2,Y2,tI);
