function [kx,kz] = thermalcracking(nx,nz,kon,koff,T)

% some constants
Apress = 0.0560;
Atemp = 7258.7;
KIc = 1e6;
Pw = 200e5;
Tve = 900;
rhoR = 3000;

% calculate cracked region
KI = Atemp*(Tve-T1)-Apress*(Pw+rhoR*g*Z);
crack = max(0,KI/KIc);

% rebuild kx and kz based on crack
kx = ones(nz,nx)*kon;
kz = ones(nz,nx)*kon;
kx(crack<1) = koff;
kz(crack<1) = koff;

