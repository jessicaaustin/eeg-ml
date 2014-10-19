% [i,d,m]=transformIDM_normalizeImages(info,data,meta) 
%
% transformIDM_normalizeTrials(info,data,meta) 
% Returns an IDM that normalizes voxel values within each image so that mean
% activity over all voxels becomes zero, and the standard deviation
% becomes 1.  Each image in the IDM is normalized independently of the others.
% More specifically, the value of each voxel v in image I gets the new value
% (v- mean(I)) / std(I).  Here mean(I) is the mean activation over all voxels in
% image I, and std(I) is the standard deviation of the voxel activities in image
% I.  In other words, this normalization resuls in an image whose mean (over all
% voxels) becomes zero, and whose standard deviation (over all voxels) becomes 1.
% Notice that the correlation between two unnormalized images I1 and I2 can be
% computed by taking the inner product of their post-normalization images.  We
% hope this normalization will be helpful in elimating effects of slow temporal
% drifts and other noise that harm Naive Bayes and other classifiers.
%
%
% Example: 
% [info2,data2,meta2] = transformIDM_normalizeImages(info,data,meta) 
%
% History
% - 7/20/2005 TMM Created file.
% - 7/25/05 TMM changed std(x,1) to std(x), after Wei convinced me to!

function [rinfo,rdata,rmeta] = transformIDM_normalizeImages(info,data,meta)
  rdata = cell(size(data));
  for trial=1:1:size(data,1)
    d=data{trial};
    rd=zeros(size(d));
    for time=1:size(d,1)
      mu=mean(d(time,:));
      dev=std(d(time,:));
      rd(time,:)= (d(time,:) - mu)/dev;
    end;
    rdata{trial}=rd;
  end
  rinfo=info;
  rmeta=meta;
 
  
