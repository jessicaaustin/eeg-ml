function diy = daysinyear(year)
%DAYSINYEAR Number of days in a year.
%
%   DAYSINYEAR(YEAR) returns the number of days in the given year.
%
%   If no year is specified, the current year is used.  Gregorian calendar is assumed.

   if nargin < 1
      clk = clock;
      year = clk(1);
   end

   diy = 365 + isleapyear(year);
