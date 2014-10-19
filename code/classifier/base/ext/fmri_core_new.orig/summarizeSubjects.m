% summarizeSubjects(study,subjects,rois,<outfile>)
% Loads data for each subject in turn, and prints the maximum and minimum
% voxel activities.  Useful for detecting suspicious data in a study.
% If the optional argument <outfile> is provided, then the output is directed
% to this file as well as the terminal screen.
% Example:
%  summarizeSubjects('data-starplus-sp',{'04772_20_sp3terC'},
%                    {'LIPS' 'LIT' 'CALC'}, 'subject-max-min-values.txt')
%
%
% History
% - 8/19/02 TMM and Xuerui - Created file.
%


function [mn,mx] = summarizeSubjects(varargin )

  l = length(varargin);
  if l < 1
    fprintf(1,'syntax: summarizeSubjects(study,subjects,rois,<outfile>)\n\n');
    fprintf(1,'- study - e.g., data-starplus-sp, data-categories,...\n');
    fprintf(1,'- subjects - cell array of subjects \n');
    fprintf(1,'- rois - cell array of ROIs to consider\n');
    fprintf(1,'- <outfile> - optional, if given then output director to this file as well as terminal screen.\n');
    fprintf(1,'\n');    
    return;
  else
    study           = varargin{1};
    subjects        = varargin{2};
    rois            = varargin{3};
    
    outfile='';
    if l > 3
      outfile = varargin{4};
    end
  end

  study
  subjects
  rois
  outfile
  
  if ~isempty(outfile) 
    fp = fopen(outfile,'w');  % write output here as well as terminal screen
  end
  
  for subj=1:1:length(subjects)
    [info,data,meta] = loadSubjectdata(study,subjects{subj},rois);
    mx=-9999;
    mn=9999;
    for j=1:1:length(data)
      mx=max(mx,max(data{j}(:)));
      mn=min(mn,min(data{j}(:)));
    end
    fprintf(1,'\n  *****  subject %s max=%f, min=%f\n',subjects{subj},mx,mn);
    if ~isempty(outfile)
      fprintf(fp,'\n  *****  subject %s max=%f, min=%f\n',subjects{subj}, ...
	      mx,mn);
    end
  end
  
  if ~isempty(outfile) 
    fclose(fp);
  end
  
  return;
  
