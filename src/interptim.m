function [z] = interptim(xi,yi,zi,x,y);
%This function replaces interp2 for simple linear interpolations
%of the thermodynamic properties computed by Will's EOS.  It only
%checks to make sure that x and y are not out of range (200 to 300
%bars for pressure and 0 to 400 C for temp)  This is a very dumb
%function -- use carefully.

%Test ranges (this is valid for hydrotab7)
%if x < 200 | x > 800; error('Pressures out of range');end;
%if y < 0 | y > 800; error('Temps out of range');end;

% deal with temps above 800 for use in uncracked region of domain
y(y>800)=800;

[nrows,ncols] = size(zi);
mx = prod(size(xi)); my = prod(size(yi));
s = 1 + (x-xi(1))/(xi(mx)-xi(1))*(ncols-1);
t = 1 + (y-yi(1))/(yi(my)-yi(1))*(nrows-1);

% Matrix element indexing
ndx = floor(t)+floor(s-1)*nrows;

% Compute intepolation parameters, check for boundary value.
d = find(s==ncols);
s(:) = (s - floor(s));
if length(d)>0, s(d) = s(d)+1; ndx(d) = ndx(d)-nrows; end

% Compute intepolation parameters, check for boundary value.
d = find(t==nrows);
t(:) = (t - floor(t));
if length(d)>0, t(d) = t(d)+1; ndx(d) = ndx(d)-1; end


% Interpolate
z =  ( zi(ndx).*(1-t) + zi(ndx+1).*t ).*(1-s) + ...
   ( zi(ndx+nrows).*(1-t) + zi(ndx+(nrows+1)).*t ).*s;
