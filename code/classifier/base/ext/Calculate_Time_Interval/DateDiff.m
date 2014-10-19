function [Days, Secs] = DateDiff(T1, T2)
% DATEDIFF Calculate elapsed days and seconds between 2 dates
% [Days, Secs] = DateDiff(T1, T2)
% INPUT:
%   T1, T2: 2 dates as accepted by Matlab's DATENUM function, e.g.:
%             - serial date numbers (e.g. from NOW)
%             - [1 x 3] and [1 x 6] date vectors (see DATEVEC)
%             - strings (see DATESTR).
%           The dates can be cell strings (Matlab >= 7.? only) or vectors. Then
%           either one is a scalar or both have the same number of elements.
%           If T2 is omitted, the current date is used.
%
% OUTPUT:
%   Days:   Number of elapsed days between the 2 dates.
%   Secs:   Number of additional seconds. 0 <= Secs < 86400.
%
%   If T1 or T1 have more than 1 element, Days and Secs are [N x 1] vectors.
%
% EXAMPLES:
%   [D, S] = DateDiff('09-Sep-1900 09:09:18', '31-Oct-2010 23:33:08')
%   % D = 40229, S = 51830
%
%   [D, S] = DateDiff([2000, 1, 1], [2010, 12, 9])
%   % D = 3995, S = 0
%
% Get the age of your family members in days:
%   D = DateDiff({'01-Jan-1970', '17-Mar-1973', '15-Jun-2001'}, now)
%
%
% See also: DATENUM, DATEVEC, ETIME.


% Initialize: ==================================================================
switch nargin
   case 1
      T2 = now;
   case 2
      % Nothing to do
   otherwise
      error(['FEX:', mfilename, ':BadNInput'], ...
         [mfilename, ': 1 or 2 inputs required']);
end

% Do the work: =================================================================
% Try to convert the dates to numbers:
try
   T1num = datenum(T1);
   T2num = datenum(T2);
catch  % No caught object to support old Matlab 6.5...
   error(['FEX:', mfilename, ':BadDate'], ...
      [mfilename, ': Cannot convert date: ', lasterr]);
end

% Calculate the positive difference -- but perhaps the sign is useful?
try
   D = abs(T2num(:) - T1num(:));
catch
   error(['FEX:', mfilename, ':SizeNotMatching'], ...
      [mfilename, ': The sizes of the input dates do not match.']);
end

Days = fix(D);
Secs = round(mod(D, 1) * 86400);

return;
