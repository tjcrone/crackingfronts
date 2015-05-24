function [Nu] = nusselt(nx,nz,d,T,P,lamdam,Tbt,Tbb,kx,kz,g, ...
   Pbt,Pbb,Pbr,Pbl)
% This function computes the nusselt number for horizontal
% layers in a 2-D porous convection computational domain that
% is heated from below.  The top and bottom boundaries should
% be maintained at a constant uniform temperature, and the box
% should be closed to fluid flow.  The sides should be insulated.
% Nusselt numbers are computed only for the boundaries that
% separate horizontal control volume bands in the model domain.

% calculate rhof,cf,qx,qz
global TT PP RHO CP BETA ALPHA
if isempty(TT)
   load('../hydrotables/hydrotab7.mat');
end

rhof = interptim(PP,TT,RHO,P./100000,T);
cf = interptim(PP,TT,CP,P./100000,T);
mu = dynvisc(T);
[qx,qz] = darcy(nx,nz,P,rhof,rhobb,kx,kz,mu,g,d,Pbt,Pbb,Pbr,Pbl,T);

% compute conductive transport
Qc = -lamdam.*(Tbb(1,1)-Tbt(1,1)).*((nx+1)/(nz+1));

% calculate the conductive heat transport across each boundary (with convection)***
Qcc = -lamdam.*(T(2:end,:)-T(1:end-1,:));

% ***Calculate the Advective Heat Transport***
Qan = T(2:end,:)  .*rhof(2:end,:)  .* cf(2:end,:)  .*qz(2:end-1,:).*(qz(2:end-1,:)<0);
Qas = T(1:end-1,:).*rhof(1:end-1,:).* cf(2:end,:)  .*qz(2:end-1,:).*(qz(2:end-1,:)>=0);
Qa = (Qan+Qas)*d;

%***Calculate Total Q***
Q = sum(Qcc+Qa,2);

%***Calculate the Nusselt Number (for each band of horizontal boxes)
Nu = Q./Qc;
