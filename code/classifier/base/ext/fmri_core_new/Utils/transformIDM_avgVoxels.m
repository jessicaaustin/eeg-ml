% transformIDM_avgVoxels(info,data,meta,<superVoxels>,<absoluteVal>)
% 
% Returns a copy of info,data, and meta, with data defined in terms of the
% specified 'superVoxels'.  The optional input argument superVoxels
% is a cell array of arrays, each of which specifies a set of
% voxels to be combined into a superVoxel. If it is not specified,
% then we combine all voxels into one superVoxel. The ith
% supervoxel is assigned a pseudo-coordinate of (i,1,1).
% The optional argument 'absoluteVal' determines whether to average the
% signed voxel values, or their absolute values.  Set it to 1 if you wish
% to average the absolute values -- it defaults to 0, which specifies
% averaging the signed values.
%
% Example:  
%  [info,data,meta] = transformIDM_avgVoxels(info,data,meta,{[3,4,5],[10,12]});
%   returns an IDM where the data contains just two voxels, one containing
%   the average of voxels 3,4 and 5, and the second the avg of voxels 10
%   and 12.
%
% Dependencies: none
%
% History
% - 9/7/02 TMM Created file.
% - 10/15/02 Xuerui now support to average all voxels into a
%                   superVoxel.
% - 10/21/02 Xuerui optimized the computation in all voxels case
% - 1/1/03 Xuerui added average active voxels(with intensity greater than 0) only
%                 right now has been commented out

function [rinfo,rdata,rmeta] = transformIDM_avgVoxels(varargin)
  l = length(varargin);
  if l < 3
    fprintf(1,'syntax: transformIDM_avgVoxels(info,data,meta,<superVoxels>,<absoluteVal>)');
    return;
  end
  info = varargin{1};
  data = varargin{2};
  meta = varargin{3};
  if l < 4
    voxelIndex = 1:1:meta.nvoxels;
    superVoxels = {voxelIndex};
  else
    superVoxels = varargin{4};
  end
  
  if l>4
    absoluteVal = varargin{5};
  else
    absoluteVal = 0;
  end			     

  rinfo= info;
  rmeta=meta;
  rmeta.nvoxels=length(superVoxels);
  rdata=cell(rmeta.ntrials,1);    

  % assign the ith supervoxel the pseudo-coordinate (i,1,1)
  for i=1:1:length(superVoxels)
    rmeta.colToCoord(i,:)=[i,1,1];
    rmeta.coordToCol(i,1,1)=i;
  end
  
  if l<4
     for i=1:1:rmeta.ntrials
%          for j = 1:1:size(data{i},1)
%              sum = 0;
%              index = find(data{i}(j,:)>0);
%              normalizer = length(index);
%              for m = index
%                  sum = sum + data{i}(j,m)/normalizer;
%              end
%              rdata{i}(j,1) = sum;
%          end
       rdata{i} = mean(data{i},2);
     end
     return;
  end
  
  % calculate the supervoxel values as avg of their component voxels
  for i=1:1:length(data) % for each trial
    rdata{i}=zeros(size(data{i},1),length(superVoxels));
    for j=1:1:size(data{i},1) % for each snapshot in the trial
      for k=1:1:length(superVoxels) % for each supervoxel in the snapshot
	normalizer=length(superVoxels{k});
	for m=1:1:length(superVoxels{k}) % for each voxel in the supervoxel 
	  if absoluteVal
	    rdata{i}(j,k)=rdata{i}(j,k) + (abs(data{i}(j,superVoxels{k}(m)))/ normalizer);
	  else
	    rdata{i}(j,k)=rdata{i}(j,k) + (data{i}(j,superVoxels{k}(m))/normalizer);
	  end
	end
      end
    end
  end
