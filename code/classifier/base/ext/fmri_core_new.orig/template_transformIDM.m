% [info,data,meta,[arg]] = template_transformIDM( info,data,meta,[arg])\n\n');
%
% Transforms an IDM (Info+Data+Meta) triple of structures into
% another one.
%
% Please consult readme.IDM in order to see further specs of
% what the transformation has to preserve or maintain consistent.
%
% In:
% - info+data+meta
% - other arguments
%
% Out:
% - info+data+meta
% - other outputs - could be information regarding the
% transformation (e.g. if you select only certain voxels of the
% IDM, you could return a list of their column numbers in the
% original IDM).

%
% History:
% - YYYY MM DD - email - what they did
%
% Notes:
% 
% Any warnings, caveats or anything else

function [transInfo,transData,transMeta,extra] = template_transformIDM( varargin )

  % varargin lets you give the function a variable number of
  % arguments, and act according to what is given. The following
  % code is an example of how to process one mandatory and one
  % optional argument (respectively arg1 and arg2)
  
  l = length(varargin);
  if l < 3
    help template_transformIDM; return;
  else
    info = varargin{1};
    data = varargin{2};
    meta = varargin{3};

    arg = 'default';
    if l > 3
      arg = varargin{4};
    end
  end
  
  % You can have one or more returns. The program that
  % called this function can pick as many returns as it wants,e.g.
  % if the call were
  % [return1] = template(arg1);
  % the first return from template.m would be used, and the second
  % ignored.
  extra = [];
  transInfo = info;
  transData = data;
  transMeta = meta;

 
% If you want, add a sequence of tests which return 0 if all went
% well or anything else if not. Eventually we might generate test
% suites for the code automatically from these functions...

function [errorStatus] = testThis()

