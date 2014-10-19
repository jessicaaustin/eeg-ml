clear all; clc;

% T = [1900 9 9 14 00 00];
% [YearRef MonthRef DayRef HourRef MinuteRef SecondRef] = [1900 9 9 14 00 00];
% YearRef = 2012; % 1900;
YearRef = 2005; % 1900;
MonthRef = 8; % 9;
DayRef = 15; % 9;
HourRef = 1; % 9;
MinuteRef = 0; % 9;
SecondRef = 0; % 18;


% futureMoment = [2012 9 9 14 00 00];

T = [YearRef MonthRef DayRef HourRef MinuteRef SecondRef];
% [year month dayOfMonth hour min sec]
[currentYear currentMonth currentDayOfMonth currentHour currentMin currentSec] = datevec(now); 
currentMoment = [2005 8 22 23 0 0];
% currentMoment = [currentYear currentMonth currentDayOfMonth currentHour currentMin currentSec];

time1 = currentMoment;
time2 = T;

earlierTime = givesTheEarlierTime (time1, time2);

if (sum(earlierTime - time1) == 0)
    laterTime = time2;
else
    laterTime = time1;
end

fprintf('The time interval between earierTime [%s] and laterTime [%s] is \n', datestr(earlierTime), datestr(laterTime));

% [totalnumberOfDays, totalnumberOfSecs] = timeDifference(time1, time2);
[totalnumberOfDays, totalnumberOfSecs] = timeDifference(time2, time1);

numberOfYears = totalnumberOfDays/365;
[hour, minute, second] = sec2hms(totalnumberOfSecs);

% determine number of years
residueYear = numberOfYears - floor(numberOfYears);
numberOfYears = floor(numberOfYears);

% determine number of days
numberOfDays = residueYear*365;
residueDay = numberOfDays - floor(numberOfDays);
numberOfDays = floor(numberOfDays);

% determine number of hours
baseOfDayInHours = 24;
numberOfHours = residueDay*baseOfDayInHours;
hour = hour + numberOfHours;

% adjust the units to allow lower unit to overflow to the upper unit, given the base
[ numberOfDays, hour ] = adjustUnits2(numberOfDays, hour, baseOfDayInHours);

if (hour >= floor(numberOfHours))
    residueHours = hour - floor(numberOfHours);
    numberOfMins = residueHours*60 + minute; 
    numberOfHours = floor(numberOfHours);
else
    numberOfHours = hour;
    numberOfMins = minute;
end

% adjust the units to allow lower unit to overflow to the upper unit, given the base
baseOfHourInMins = 60;
[ numberOfHours, numberOfMins ] = adjustUnits2(numberOfHours, numberOfMins, baseOfHourInMins);

% determine number of mins
residueMins = numberOfMins - floor(numberOfMins);
numberOfSecs = residueMins*60 + second; 
numberOfMins = floor(numberOfMins);

% determine number of secs
% adjust the units to allow lower unit to overflow to the upper unit, given the base
baseOfMinInSecs = 60;
[ numberOfMins, numberOfSecs ] = adjustUnits2(numberOfMins, numberOfSecs, baseOfMinInSecs);

fprintf('%s years, %s days, %s hours, %s mins, %s secs.\n', ...
    num2str(numberOfYears), num2str(numberOfDays), num2str(numberOfHours), num2str(numberOfMins), num2str(numberOfSecs));

% Jan Simon's rendition
[D, S] = DateDiff('15-Aug-2005 01:00:00', '22-Aug-2005 23:00:00');
fprintf('\n%s days, %s secs.\n', ...
    num2str(D), num2str(S));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
Notes:
Calculate Time Interval
--------------------------------------
main executing reference usage: usage_timeDifference.m 

Sample output:
The time interval between earierTime [09-Sep-1900 09:09:18] and laterTime [31-Oct-2010 23:33:08] is 
110 years, 79 days, 14 hours, 23 mins, 50.17 secs.


The objective is to compute the time interval taking into account the leap days that are subsumed within 
the time interval. 

Note:
    The logic subsumes the time interval by 
        yearEnd, monthEnd, dayEnd - timeEarlier (within that year)
        + time interval (excluding that 2 years)
        + timeLater - yearEnd, monthEnd, dayEnd (within that year)

It differentiates the earlier time and later time, ie. timeDifference (earlierTime, laterTime) gives the 
same result as timeDifference (laterTime, earlierTime).

leapDayHasPassed.m illustrates a possible approach to determine if a certain day of concern is already over 
at a given time.

Caveat: does not include leap secs 
%}

%{
% currentMoment = futureMoment;
currentYear = 2012;
currentMonth =  9;
currentDayOfMonth =  19;
currentHour = 14;
currentMin = 0;
currentSec = 0;
%}