function [z] = interptim(xi,yi,zi,x,y)
% This function replaces interp2 for faster linear interpolations.
% This is a very dumb function so use carefully.

% deal with temps above 800 for use in uncracked region of domain
y(y>800)=800;

% setup
[nrows,ncols] = size(zi);
mx = prod(size(xi)); my = prod(size(yi));
s = 1 + (x-xi(1))/(xi(mx)-xi(1))*(ncols-1);
t = 1 + (y-yi(1))/(yi(my)-yi(1))*(nrows-1);

% matrix element indexing
ndx = floor(t)+floor(s-1)*nrows;

% compute intepolation parameters, check for boundary value
d = find(s==ncols);
s(:) = (s - floor(s));
if length(d)>0, s(d) = s(d)+1; ndx(d) = ndx(d)-nrows; end

% compute intepolation parameters, check for boundary value
d = find(t==nrows);
t(:) = (t - floor(t));
if length(d)>0, t(d) = t(d)+1; ndx(d) = ndx(d)-1; end

% interpolate
z =  ( zi(ndx).*(1-t) + zi(ndx+1).*t ).*(1-s) + ...
   ( zi(ndx+nrows).*(1-t) + zi(ndx+(nrows+1)).*t ).*s;
