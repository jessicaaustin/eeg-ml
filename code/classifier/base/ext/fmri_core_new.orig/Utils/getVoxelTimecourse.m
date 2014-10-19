% getVoxelTimecourse(info,data,meta,x,y,z)
%
%  Returns a column vector containing the entire time course of voxel
%  throughout the IDM
%
% WISH LIST: add error checking for whether there is data for x,y,z
%  
% History: 
% 8/13/02 Tom - created.  

function [vtc] = getVoxelTimecourse(info,data,meta,x,y,z)
  tcLength = sum([info.len]);
  vtc=zeros(tcLength,1);
  col = meta.coordToCol(x,y,z);
  beginPos=1;
  for i=1:1:meta.ntrials
    trial=data{i};
    voxTrial=trial(:,col);
    endPos=beginPos-1+info(i).len;
    vtc(beginPos:endPos)=voxTrial;
    beginPos=endPos+1;
  end
