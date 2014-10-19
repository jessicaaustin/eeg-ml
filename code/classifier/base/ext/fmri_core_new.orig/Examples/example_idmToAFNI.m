% Example of how to use idmToAFNI
%
% WARNING: run this in a clean directory, as it generates several files
%

subject = '05886';
study   = 'data-categories';
information = localInformation('data-categories');
%rois    = information.rois;
rois    = {'LIT' 'CALC'};

% load the data - all ROIs merged into a single one
[info,data,meta] = loadSubjectdata(study,subject,rois);


prefix  = 'prefix'; % prefix of files that are output
trial   = 2; % trial that is to be plotted
img     = 1; % image # within trial (if 0, do a 4D plot)
example = 3; % which demo to run (see comments in each case statement)

switch example
  
 case {1}
  % vanilla version  
  idmToAFNIfiles_wrapper(info,data,meta,subject,study,prefix,trial,img);
  
 case {2}
  % Scale data
  % By default, data is not scaled at all (see vanilla version).
  % You can suggest a scale by supplying two additional parameters,
  
  % Scale by given range (e.g. [-3 3]) - used by code that needs to
  % compare several conditions, for instance
%  idmToAFNIfiles_wrapper(info,data,meta,subject,study,prefix,trial,img,'scale',[-3,3]);
  
  % Scale by the 1% and 99% percentiles in the data
%  idmToAFNIfiles_wrapper(info,data,meta,subject,study,prefix,trial,img,'scale','percentiles');
  
  % Scale by the min and max in the data - very volatile
%  idmToAFNIfiles_wrapper(info,data,meta,subject,study,prefix,trial,img,'scale','minmax');

 case {3}
  % Tag data
  % You can create tags to mark the images. The first element in
  % the labels cell array will appear in the slice with the largest
  % z. The first element in a cell appears on top, the second one
  % at the bottom (less space)
  %
  % You may want to try assembling a multi-slice view. Do it by
  % pressing "Mont" on a axial window, then picking 4 across, 4
  % down (we have 4x4=16 slices) and spacing = 4 (there's 4
  % anatomical slices for each functional one). Then adjust the
  % slider until all the functional slices are in view.
  
  labels{1} = {subject trial};
  labels{2} = 'twinkle';
  labels{3} = 'twinkle';
  labels{4} = 'little';
  labels{5} = 'star';
  idmToAFNIfiles_wrapper(info,data,meta,subject,study,prefix,trial,img,'labels',labels);
end