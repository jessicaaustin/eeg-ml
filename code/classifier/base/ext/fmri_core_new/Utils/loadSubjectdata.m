%loadSubjectdata(Study,Subject,ROIs)
% loads data for this study and subject, for the specified rois, from the
% appropriate files.
% Study is the name of the fMRI study (e.g., 'data-sentence-ps')
% Subject is ID of subject (e.g., '055868')
% ROIs is a cell array of ROI names (e.g., {'LIT' 'CALC'})
% Returns the Info,Data, and Meta datastructures that represent the dataset.
% Example:  
%  [info,data,meta] = loadSubjectdata('data-categories', '05581', {'LIT' 'CALC'})
%
% Dependencies:
% - process_dataToExamples depends on
%
% History
% - 020404 - fp
% - Added support for the categories data set and made example code
% standalone (rather than functions callable inside this)
% - Aug13,2002 Tom - added fields to meta: subject, study, roi
% - Aug18,2002 Tom - made separateROIs an optional argument
% - Aug24,2002 fp  - a cell array of ROIs is returned, even if we
%                    ask for a single one -> allows the rest of the
%                    code to be the same for any number of ROIs.
% - Aug26,2002 tm - undid the Aug 24 change.  loadSubjectdataMult provides
%                   that functionality.  See the master documentation on
%                   that function.


function [rinfo,rdata,rmeta] = loadSubjectdata( varargin )
  l = length(varargin);
  
  if l < 1
    fprintf(1,'syntax: loadSubjectdata(study,subject,rois)\n\n');
    fprintf(1,'- study - data-starplus-sp/ps or data-categories\n');
    fprintf(1,'- subject - which subject inside the study\n');
    fprintf(1,'- rois - which ROIs to use\n');
    fprintf(1,'\n');    
    return;
  else
    study           = varargin{1};
    subject         = varargin{2};
    rois            = varargin{3};
  end

  % let the work be done by loadSubjectdataMult
  [is,ds,ms]=loadSubjectdataMult(study,subject,rois,0);
  rinfo=is{1};
  rdata=ds{1};
  rmeta=ms{1};
