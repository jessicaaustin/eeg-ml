% Reads in voxel coordinates that define the ROIs, creating two new fields
% for the 'meta' datastructure:
%  meta.rois: [1x nrois] array of structures, each containing subfields:
%     name: the name of the ROI, e.g., 'LIT'
%     coords: [n x 3] array with the coordinates of all n voxels in this ROI
%     columns: [1 x n]  array with the column numbers of all n voxels in this ROI
%  meta.colToROI: a {nvoxels x 1} cell array containing the name of the
%     ROI for each voxel.  That is, meta.colToROI{n} contains the name of the
%     ROI for the voxel in column n.  e.g., meta.colToROI{20}='lit'
%
% Example:
% [i,d,m]=createColToROI(i,d,m);
%
% Efficiency: perhaps it would be more efficient if meta.colToROI were a
% normal array containing integer id's of the ROIs instead of their names
% stored as strings...
%
% Depends on: loadROIcoords.m
%
% History
% 9/13/02 Tom - Created.
%
% ROIrootdir='/shared/fmri/data/data-brainlex/08057/mroi/';



function [i,d,m] = createColToROI(i,d,m)

  % if m.rois not present, then need to load ROIs first
  if ~isfield(m,'rois')
      fprintf('\nloading ROI coordinates to create meta.rois and meta.colToROI\n');
      [i,d,m]=loadROIcoords(i,d,m);
  end

  % now create m.colToROI, and create m.rois(j).columns
  m.colToROI=cell(m.nvoxels,1);
  for j=1:1:length(m.rois) % for each ROI
    ROIname=m.rois(j).name;
    m.rois(j).columns=[];  
    for k=1:1:size(m.rois(j).coords,1); % each coord of that ROI
      xyz=m.rois(j).coords(k,:);
      %% IMPORTANT!!! Francisco says that each voxel coordinate in the
      % original data begins with 0, and he adds 1 to each coordinate so
      % that Matlab can instead have them all begin with the value 1.
      % Hence, we do the same here:
      col=m.coordToCol(xyz(1)+1,xyz(2)+1,xyz(3)+1);
      if col>0
	m.colToROI{col}=ROIname;
	m.rois(j).columns= [m.rois(j).columns, col];
	%	  col
	%	  xyz
      end
    end
  end
  
  
