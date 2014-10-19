% load ROI coordinate information from subject data directory, adding it
% to the meta of the info/data/meta.
%
% creates meta.rois as a cell array where each cell has the properties
% 'name' and 'coords'.  'name' is the string name of the subdirectory
% corresponding to the ROI.  'coords' is a [numVoxels x 3] array giving
% the x,y,z coordinate for each voxel in the ROI.  
%
% Example:
% [i,d,m]=loadROIcoords(i,d,m,'/shared/fmri/data/data-brainlex/08057/mroi/');
%
% History
% 9/13/02 Tom - Created.
% 11/21/02 Tom - made directory paths work under Windows and Linux
% 3/15/04 Indra - modified to consider only rois available through mri_reference

function [i,d,m] = loadROIcoords(i,d,m)
  
  information = localInformation(m.study);
  rois = information.rois;
  ROIrootdir=fullfile(information.dataplace, m.study, m.subject,'mroi');
  %dirs=dir(ROIrootdir);
  roinum=0;
  for j=1:1:length(rois)
    roi = rois(j);
    roinum=roinum+1;
    m.rois(roinum).name=roi{1};
    ROIdir=fullfile(ROIrootdir,roi{1},'allvoxels');
    m.rois(roinum).coords=load(ROIdir, '-ascii');
    end
  end
  
  
