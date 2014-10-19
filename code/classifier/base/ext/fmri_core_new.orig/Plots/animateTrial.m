% animateTrial(i,d,m,trialNum,z)
%
% Create a movie of brain activity for one z slice of trial trialNum of the
% IDM, and display it.  Returns the movie matrix.
%
% Example: Show and return a movie for brain slice z, trial trialNum
%  M=animateTrial(i,d,m,trialNum,z);  
%  movie(M,3,2);     % replay it 3 times, at 2 frames/second
% 
% History
% - Aug17,2002 Tom - created.


function [M] = animateTrial(i,d,m,trialNum,z)
  img=zeros(m.dimx,m.dimy);
  
  minActiv=-10.0;
  maxActiv=10.0;
  
  for j=1:1:i(trialNum).len
    % for each snapshot in the trial, create and store a frame
    img=snapToImage(i,d,m,trialNum,z,j,minActiv,maxActiv);
    image(img);
    axis image;
    M(j) = getframe;
  end
  
  % finally, play the movie at 1 frame per sec (assume TR=1sec)
  movie(M,1,2);
  
