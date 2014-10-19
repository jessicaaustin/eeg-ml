% animate6Trial(i,d,m,trialNum,z)
%
% Create a 4-slice movie of brain activity for trial trialNum of the IDM, and
% display it.  Displays slices z,z+1,z+2,z+3. Returns the movie matrix.
%
% Example
%  M=animate6Trial(i,d,m,trialNum,z);  
%  movie(M,2,3);   % replay it twice, as 3 frames per second.
% 
% History
% - Aug17,2002 Tom - created.


function [M] = animate6Trial(i,d,m,trialNum,z)
  img=zeros(m.dimx,m.dimy);
  
  minActiv=-10.0;
  maxActiv=10.0;
  
  bigimg=zeros(128,128);
  for j=1:1:i(trialNum).len
    % for each snapshot in the trial, create and store a frame for
    % each z slice
    bigimg(1:64,1:64)=snapToImage(i,d,m,trialNum,z,j,minActiv,maxActiv);
    bigimg(1:64,65:128)=snapToImage(i,d,m,trialNum,z+1,j,minActiv,maxActiv);
    bigimg(65:128,1:64)=snapToImage(i,d,m,trialNum,z+2,j,minActiv,maxActiv);
    bigimg(65:128,65:128)=snapToImage(i,d,m,trialNum,z+3,j,minActiv,maxActiv);
    bigimg(129:192,1:64)=snapToImage(i,d,m,trialNum,z+4,j,minActiv,maxActiv);
    bigimg(129:192,65:128)=snapToImage(i,d,m,trialNum,z+5,j,minActiv,maxActiv);
    
    image(bigimg)
    axis image;
    M(j) = getframe;
  end
  
  % finally, play the movie at 1 frame per sec 
  movie(M,1,2);
  
