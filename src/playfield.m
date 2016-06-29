function [] = playfield(inputfield, pauselength)
% This function plays the field that is passed, usually Tout.
%
% Timothy Crone (tjcrone@gmail.com)


% set up figure
figure;
im1 = imagesc(inputfield(:,:,1));
set(gca,'dataaspectratio',[1 1 1]);
colormap(jet);
colorbar;

% get length of input field that is filled (non-zero)
inputlength = max(find(squeeze(max(max(inputfield)))));

% get range of input field

% loop
for i=1:inputlength
  set(im1,'cdata',inputfield(:,:,i))
  title(i)
  drawnow
  caxis([0 600])
  pause(pauselength)
end
