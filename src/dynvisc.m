function [mu] = dynvisc(T)
%This function computes the dynamic viscosity of pure water as
%a function of temperature. Viscosity is returned in units of
%Pascal seconds. T can be a scalar, vector or matrix of any
%size or shape. The valid range of temperatures for this
%function is 0-400 degrees centigrade. Pressure for the
%temperature range 0-100 is approximately 1 bar. Equations
%are taken from the following reference:
%
%S. Schoofs, U.Hansen, Depletion of a brine layer at the base of
%ridge-crest hydrothermal systems, Earth Planet. Sci. Lett. 180
%(2000) 341-353

%Error Checking
if max(max(max(T)))>1000 || min(min(min(T)))<0
   error('Temperature out of valid range for viscosity function.');
end

%Initialize Output
mu = zeros(size(T));

%Find Locations >100, <100
great100 = find(T>100);
less100 = find(T<=100);

%Compute Viscosity
mu(less100) = ((((T(less100)-20)*1.551e-2)+1).^(-1.572))*1e-3;
mu(great100) = (2.414e-5)*10.^(247.8./(T(great100)+133.15));
