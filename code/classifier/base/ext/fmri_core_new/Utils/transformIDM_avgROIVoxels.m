% transformIDM_avgROIVoxels(info,data,meta,<superVoxels>,<absoluteVal>)
% 
% transformIDM_avgROIVoxels(info,data,meta,<ROIs>)
% Returns a copy of info,data, and meta, with each ROI replaced by a single
% 'superVoxel', whose activation is the mean activation of voxels in that ROI.
% ROIs is a cell array specifying the text names of the ROIs to include.  If ROIs
% is not given, then all ROIs in meta.rois are used.  The ith ROI supervoxel is assigned a
% pseudo-coordinate of x=i,y=1,z=1.  (implemented using transformIDM_avgVoxels).
% **NOTE: this function assumes you have run createColToROI beforehand ***
% Example:  
%  [i,d,m] = transformIDM_avgROIVoxels(i,d,m,{'LT' 'RB'});
%   returns an IDM where data contains just two voxels, the first
%   giving the spatial average of voxels in LT, the second giving the mean activity of voxels in RB.
%
% Dependencies: uses  transformIDM_avgVoxels
%
% History 
% - 2/03 TMM Created file.
% -10/04 TMM  added code to add meta.colToROIsupervoxel=ROIs;  This
%                      adds information to meta to show which
%                      column corresponds to which ROI.


function [rinfo,rdata,rmeta] = transformIDM_avgROIVoxels(varargin)
  l = length(varargin);
  if l < 3
    fprintf(1,'syntax: transformIDM_avgVoxels(info,data,meta,<ROIs>)');
    return;
  end
  info = varargin{1};
  data = varargin{2};
  meta = varargin{3};
  if l < 4
   ROIs={meta.rois.name};
  else
   ROIs = varargin{4};
  end     
  
 % create the supervoxels cell array from the ROI names
  supervoxels={};

  if ~isfield(meta,'rois')
      [info,data,meta]=createColToROI(info,data,meta);
  end
  for i=1:length(ROIs)
    roi=getROI(meta,ROIs(i));
%    i
    if (~isstruct(roi) | length(roi.columns)==0)
      fprintf('warning in transformIDM_avgROIVoxels: cannot find ROI %s  in meta', ROIs{i});
      break;
    else
      supervoxels{length(supervoxels)+1}= roi.columns;
    end
  end

  % finally create the averaged I,D,M
  [rinfo,rdata,rmeta]= transformIDM_avgVoxels(info,data,meta, ...
                                              supervoxels);
  rmeta.colToROIsupervoxel=ROIs;
  

  
