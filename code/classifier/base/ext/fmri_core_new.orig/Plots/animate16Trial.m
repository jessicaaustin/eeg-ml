% animate16Trial(i,d,m,trialNum)
%
% Create a movie of brain activity for trial trialNum of the IDM,
% and display it.  Returns the movie matrix.  Z slices are laid out
% in format:
% 1  2  3  4
% 5  6  7  8
% 9 10 11 12
%13 14 15 16
%
% Example
%  M=animate16Trial(i,d,m,trialNum);  
%  movie(M,3,2);   % replay it 3 times at 2 frames/second
% 
% History
% - Aug17,2002 Tom - created.


function [M] = animate16Trial(i,d,m,trialNum)
  img=zeros(m.dimx,m.dimy);
  
  minActiv=-10.0;
  maxActiv=10.0;
  
  bigimg=zeros(256,256);
  image(bigimg)
  axis image;
  M(1)=getframe;
  for j=1:1:i(trialNum).len
    % for each snapshot in the trial, create and store a frame for
    % each z slice
    bigimg(1:64,1:64)=snapToImage(i,d,m,trialNum,1,j,minActiv,maxActiv);
    bigimg(1:64,65:128)=snapToImage(i,d,m,trialNum,2,j,minActiv,maxActiv);
    bigimg(1:64,129:192)=snapToImage(i,d,m,trialNum,3,j,minActiv,maxActiv);
    bigimg(1:64,193:256)=snapToImage(i,d,m,trialNum,4,j,minActiv,maxActiv);
    bigimg(65:128,1:64)=snapToImage(i,d,m,trialNum,5,j,minActiv,maxActiv);
    bigimg(65:128,65:128)=snapToImage(i,d,m,trialNum,6,j,minActiv,maxActiv);
    bigimg(65:128,129:192)=snapToImage(i,d,m,trialNum,7,j,minActiv,maxActiv);
    bigimg(65:128,193:256)=snapToImage(i,d,m,trialNum,8,j,minActiv,maxActiv);
    bigimg(129:192,1:64)=snapToImage(i,d,m,trialNum,9,j,minActiv,maxActiv);
    bigimg(129:192,65:128)=snapToImage(i,d,m,trialNum,10,j,minActiv,maxActiv);
    bigimg(129:192,129:192)=snapToImage(i,d,m,trialNum,11,j,minActiv,maxActiv);
    bigimg(129:192,193:256)=snapToImage(i,d,m,trialNum,12,j,minActiv,maxActiv);
    bigimg(193:256,1:64)=snapToImage(i,d,m,trialNum,13,j,minActiv,maxActiv);
    bigimg(193:256,65:128)=snapToImage(i,d,m,trialNum,14,j,minActiv,maxActiv);
    bigimg(193:256,129:192)=snapToImage(i,d,m,trialNum,15,j,minActiv,maxActiv);
    bigimg(193:256,193:256)=snapToImage(i,d,m,trialNum,16,j,minActiv,maxActiv);
    
    image(bigimg)
    axis image;
    M(j+1) = getframe;
  end
  
  
  % finally, play the movie at 2 frames per sec
  movie(M,1,2);
  
