% snapToImage(i,d,m,trialNum,z,n,minActiv,maxActiv)
%
% Return an image from the idm data, and display it in Matlab, for
% slice z of the n-th snapshot of trial number trialNum.  If you
% want to pass in the maximum and minimum activities (i.e., the
% activity levels that will display as intensities 1 and 128), then
% set maxActiv, minActiv accordingly.  If you set these both to
% zero, they will get computed from the data.
%
% Note Matlab images require intensities in range [1,128]
% To display this in Matlab, do this:
%  image(snapToImage(i,d,m,trialNum,z,n,0,0));  axis image;
%  colormap('default');
% 
% Note Matlab function imread will read in an image (e.g., structural
% image, in case we want to do this in the future...)
%
% History
% - Aug17,2002 Tom - created.


function [img] = snapToImage(i,d,m,trialNum,z,n,minActiv,maxActiv)
  img=zeros(m.dimx,m.dimy);
  
  % get the snapshot
  snap = d{trialNum}(n,:);
  
  if maxActiv==minActiv
    maxActiv=max(snap); minActiv=min(snap);
  end
  
%  scale values to interva%l [1,128]
%  snap = snap - minActiv;
%  snap = (snap * 127/(maxActiv-minActiv)) + 1;

% scale values to interval [1,64], to use default matlab colormap
  snap = snap - minActiv;
  snap = (snap * 63/(maxActiv-minActiv)) + 1;

  % for columns (voxels) that have right z value, add to img
  for j=1:1:length(snap)
    coord=m.colToCoord(j,:);
    if coord(3) == z
      img(coord(1),coord(2))=snap(j);      
    end
  end
  
  % add in a color bar at x=0, y=[1:128]
  for j=1:1:64
    img(1:2,j)=j;
  end
  
  % these are to sanity check the color bar and axes
  img(62,:)=127;
  img(63,:)=1;
  img(64,:)=64;
  img(64,64)=1;
  img(63,64)=127;
  img(1,1)=127;

  img=img';
  
  % do this to display the image
  % image(img);
  % axis image;
  
  
