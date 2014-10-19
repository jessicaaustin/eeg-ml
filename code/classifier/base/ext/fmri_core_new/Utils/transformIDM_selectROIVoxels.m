% transformIDM_selectROIVoxels(info,data,meta,ROIs)
% 
% Returns a copy of info,data, and meta, selecting only voxels that belong to
% ROIs.  Input ROIs is a cell array of strings naming the ROIs.
%
% Example:  
%  [info,data,meta] = transformIDM_selectROIVoxels(info,data,meta,{'CALC' 'LIT'})
%   returns an IDM where the data contains just the voxels belonging to
%   CALC and LIT
%
% Dependencies: uses transformIDM_avgVoxels
%
% History
% - 2/19/03 TMM Created file.
% - 12/9/03 indra Correction in the warning printf

function [rinfo,rdata,rmeta] = transformIDM_selectROIVoxels(info,data,meta,ROIs)
  
 % create the supervoxels cell array from the ROI names
  voxelcolumns=[];
  if ~isfield(meta,'rois')
      [info,data,meta]=createColToROI(info,data,meta);
  end
  for i=1:length(ROIs)
    roi=getROI(meta,ROIs(i));
    if (~isstruct(roi) | length(roi.columns)==0)
      fprintf('warning in transformIDM_selectROIVoxels: cannot find ROI %s  in meta', ROIs{i});
      break;
    else
      voxelcolumns=[voxelcolumns, roi.columns];
    end
  end

  % finally create the I,D,M with the selected voxels
  [rinfo,rdata,rmeta]= transformIDM_selectVoxelSubset(info,data,meta,voxelcolumns);

  
