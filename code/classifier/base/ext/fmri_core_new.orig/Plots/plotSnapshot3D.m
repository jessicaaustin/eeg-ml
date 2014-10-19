% plotSnapshot3D(i,d,m,trialNum,n,minActiv,maxActiv)
%
% Same as plotSnapshot, but displays up to 16 z slice images at once.
% 
% Create a 16 slice set of images of brain activity for trial trialNum,
% snapshot n, of the IDM.  minActiv and maxActiv determine the minimum and
% maximum activation levels displayed in the color scale of the image.  If
% minActiv = maxActiv, then both will be calculated from the image.
% Displays the image, and returns it.
%
% Z slices are laid out in format: 
% 1  2  3  4 
% 5  6  7  8 
% 9 10 11 12
%13 14 15 16
%
% Example: display first snapshot of second trial.  Let program pick min
% and max activations.
%  img = plotSnapshot3D(i,d,m,2,1,0,0);  
% 
% History
% - 9/5/02 Tom - created.


function [bigimg] = plotSnapshot3D(i,d,m,trialNum,n,minActiv,maxActiv)
  img=zeros(m.dimx,m.dimy);
  
  % get the snapshot
  snap = d{trialNum}(n,:);
  
  if maxActiv==minActiv
    maxActiv=max(snap); minActiv=min(snap);
  end
  
  bigimg=zeros(256,256);
  image(bigimg)
  axis image;

  % for each z slice, create the image and add it to the bigimg array
  bigimg(1:64,1:64)=snapToImage(i,d,m,trialNum,1,n,minActiv,maxActiv);
  bigimg(1:64,65:128)=snapToImage(i,d,m,trialNum,2,n,minActiv,maxActiv);
  bigimg(1:64,129:192)=snapToImage(i,d,m,trialNum,3,n,minActiv,maxActiv);
  bigimg(1:64,193:256)=snapToImage(i,d,m,trialNum,4,n,minActiv,maxActiv);
  bigimg(65:128,1:64)=snapToImage(i,d,m,trialNum,5,n,minActiv,maxActiv);
  bigimg(65:128,65:128)=snapToImage(i,d,m,trialNum,6,n,minActiv,maxActiv);
  bigimg(65:128,129:192)=snapToImage(i,d,m,trialNum,7,n,minActiv,maxActiv);
  bigimg(65:128,193:256)=snapToImage(i,d,m,trialNum,8,n,minActiv,maxActiv);
  bigimg(129:192,1:64)=snapToImage(i,d,m,trialNum,9,n,minActiv,maxActiv);
  bigimg(129:192,65:128)=snapToImage(i,d,m,trialNum,10,n,minActiv,maxActiv);
  bigimg(129:192,129:192)=snapToImage(i,d,m,trialNum,11,n,minActiv,maxActiv);
  bigimg(129:192,193:256)=snapToImage(i,d,m,trialNum,12,n,minActiv,maxActiv);
  bigimg(193:256,1:64)=snapToImage(i,d,m,trialNum,13,n,minActiv,maxActiv);
  bigimg(193:256,65:128)=snapToImage(i,d,m,trialNum,14,n,minActiv,maxActiv);
  bigimg(193:256,129:192)=snapToImage(i,d,m,trialNum,15,n,minActiv,maxActiv);
  bigimg(193:256,193:256)=snapToImage(i,d,m,trialNum,16,n,minActiv,maxActiv);

  tstr=sprintf('%s subject%s, trial number %d, snapshot number %d', m.study,m.subject,trialNum,n);
  tstr=sprintf('%s\nregion %s\n',tstr,m.roi);
  tstr=sprintf('%s, condition=%d ', tstr,i(trialNum).cond);
  tstr=sprintf('%s amplitude=[%1.1f %1.1f]', tstr, minActiv,maxActiv);
  
  image(bigimg)
  axis image;
  title(tstr);

