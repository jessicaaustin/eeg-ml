% getROI(meta,roiNameString)
% 
% Return the struct from meta.rois that corresponds to 'roiNameString'
% Example: getROI(meta,'LT')
%
% History:  created 2/19/03, TMM
%

function [roi] = getROI(meta, roiNameStr)
  roi=0;
  for i=1:length(meta.rois)
    if strcmpi(roiNameStr,meta.rois(i).name)
      roi=meta.rois(i);
      return
    end
  end
  
      
  
