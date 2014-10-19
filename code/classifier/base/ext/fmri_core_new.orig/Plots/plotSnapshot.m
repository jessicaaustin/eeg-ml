% plotSnapshot(i,d,m,trialNum,z,n,minActiv,maxActiv)
%
% Plot an image showing slice z of the n-th snapshot of trial number
% trialNum.  If you want to pass in the maximum and minimum activities (i.e.,
% the activity levels that will display as intensities 1 and 128), then set
% maxActiv, minActiv accordingly.  If you set these both to zero, they will
% get computed from the data.
%
% Example: display slice z=3 of snapshot 4 of trial 2
%  plotSnapshot(i,d,m,2,3,4,0,0)
%
% Note Matlab images require intensities in range [1,128]
% 
% Note Matlab function imread will read in an image (e.g., structural
% image, in case we want to do this in the future...)
%
% History
% - Aug18,2002 Tom - created.


function [img] = plotSnapshot(i,d,m,trialNum,z,n,minActiv,maxActiv)
  img=snapToImage(i,d,m,trialNum,z,n,minActiv,maxActiv);
  % do this to display the image
  clf;
  colormap('default');
  image(img);
  axis image;
  
  
