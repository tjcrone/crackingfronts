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

% set up figure
figure;
im1 = imagesc(inputfield(:,:,1));
set(gca,'dataaspectratio',[1 1 1]);
colormap(jet);
colorbar;

% get length of input field that is filled (non-zero)
inputlength = max(find(squeeze(max(max(inputfield)))));

% get range of input field
fmax = double(max(max(max(inputfield))));
fmin = double(min(min(min(inputfield))));

% loop
for i=1:interval:inputlength
  set(im1,'cdata',inputfield(:,:,i))
  title(i)
  drawnow
  caxis([fmin fmax])
  pause(pauselength)
end
