function [] = playfield(inputfield, pauselength, interval)
% This function plays the field that is passed, usually Tout.
%
% Timothy Crone (tjcrone@gmail.com)

% set default pauselength
if ~exist('pauselength', 'var')
  pauselength = 0;
end
if isempty(pauselength)
  pauselength = 0;
end

% set default interval
if ~exist('interval', 'var')
  interval = 1;
end
if isempty(interval)
  interval = 1;
end

% get length of input field that is filled (non-zero)
inputlength = max(find(squeeze(max(max(inputfield)))));

% get range of input field
fmax = double(max(max(max(inputfield))));
fmin = double(min(min(min(inputfield))));

% set up figure
figure;
im1 = imagesc(inputfield(:,:,1));
set(gca,'dataaspectratio',[1 1 1]);
caxis([fmin fmax])
colormap(jet);
colorbar;
title('1');
drawnow;
pause(pauselength);

% loop
for i=1+interval:interval:inputlength
  set(im1,'cdata',inputfield(:,:,i));
  title(i);
  drawnow;
  pause(pauselength);
end
