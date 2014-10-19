function [totalnumberOfDays, totalnumberOfSecs] = timeDifference(time1, time2)
%TIMEDIFFERENCE 
currentMoment = time1;
T = time2;

numberOfDaysWithoutLeapYears = [31 28 31 30 31 30 31 30 30 31 30 31];

% Caveat: TimeRef and currentTime are to be determined of earlier and
% later, currency of time is not relevant, thus is a misnomer
YearRef = T(1);
MonthRef = T(2);
DayRef = T(3);
HourRef = T(4);
MinuteRef = T(5);
SecondRef = T(6);

currentYear = currentMoment(1);
currentMonth = currentMoment(2);
currentDayOfMonth = currentMoment(3);
currentHour = currentMoment(4);
currentMin = currentMoment(5);
currentSec = currentMoment(6);

% check if T is past or present
elapsedTime = etime(currentMoment, T); % returns seconds

if (elapsedTime > 0)        % isPast
    TimeRef = 'past';
elseif (elapsedTime == 0)   % isPresent
    TimeRef = 'present';
else
    TimeRef = 'future';     % isFuture
end

% get the leap years within the interval
switch (TimeRef)
    case 'past'
        count = 0;
        for i = YearRef:currentYear
            count = count +1;
            leapYears(count) = isleapyear(i);
        end
        % determine from the 2 point of time has the leap day passed?
        % [yes/no? yes/no?]        
        % account for the leap day for the 2 years        
        if (isleapyear(YearRef))
            leapYears(1) = ~leapDayHasPassed( MonthRef, dayRef);   % complement because days have to be added with ref to the year end
        end
        
        if (isleapyear(currentYear))
            leapYears(end) = leapDayHasPassed( currentMonth, currentDayOfMonth); % non-complement because days have to be added with ref to the year start
        end
    case 'present'        
            leapYears = 0; % not relevant to account for leap days, since it is 'now'
    case 'future'
        count = 0;
        for i = currentYear:YearRef
            count = count +1;
            leapYears(count) = isleapyear(i);
        end
        % determine from the 2 point of time has the leap day passed?
        % [yes/no? yes/no?]        
        % account for the leap day for the 2 years
        
        if (isleapyear(currentYear))
            leapYears(1) = ~leapDayHasPassed( currentMonth, currentDayOfMonth);  % complement because days have to be added with ref to the year end
        end
        
        if (isleapyear(YearRef))
            leapYears(end) = leapDayHasPassed( MonthRef, DayRef); % non-complement because days have to be added with ref to the year start
        end
end

% getTheNumberOfDaysBetweenThat2Years, not including days within the 2 year themselves
switch (TimeRef)
    case 'past'
        % YearRef:currentYear        
        numberOfDaysBetween2YearsInterval = sum(leapYears(2:(end-2))) + (365*(abs(currentYear-YearRef)-1));
        
        numberOfDaysInThatYearToYearEnd = leapYears(1) + ...
                sum(numberOfDaysWithoutLeapYears((MonthRef+1):12)) + ...
                (numberOfDaysWithoutLeapYears(MonthRef)- DayRef);
        numberOfDaysInThisYearFromYearStart = leapYears(end) + ...
                sum(numberOfDaysWithoutLeapYears(1:(currentMonth-1))) + ...
                currentDayOfMonth;
            
        totalnumberOfDays = numberOfDaysInThatYearToYearEnd + ...
                            numberOfDaysBetween2YearsInterval + ...
                            numberOfDaysInThisYearFromYearStart;
        
    case 'present'
        totalnumberOfDays = 0;
        
    case 'future'
        % currentYear:YearRef
        numberOfDaysBetween2YearsInterval = sum(leapYears(2:(end-2))) + (365*(abs(currentYear-YearRef)-1));

        numberOfDaysInThisYearToYearEnd = leapYears(1) + ...
                sum(numberOfDaysWithoutLeapYears((currentMonth+1):12)) + ...
                (numberOfDaysWithoutLeapYears(currentMonth)- currentDayOfMonth);
        numberOfDaysInThatYearFromYearStart = leapYears(end) + ...
                sum(numberOfDaysWithoutLeapYears(1:(MonthRef-1))) + ...
                DayRef;
            
        totalnumberOfDays = numberOfDaysInThisYearToYearEnd + ...
                            numberOfDaysBetween2YearsInterval + ...
                            numberOfDaysInThatYearFromYearStart;            
end

% numberOfYears = totalnumberOfDays/365;

totalNumberOfSecsIn1Day = 1440*60; % 1440 mins

% getTheSecsOfTheParticularDaysOfThat2Years. * refer notes for common principle
switch (TimeRef)
    case 'past'
        % YearRef:currentYear    
        numberOfSecsOfDayInEarlierTime = totalNumberOfSecsIn1Day - (HourRef*60*60) - (MinuteRef*60) - SecondRef;
        numberOfSecsOfDayInLaterTime = (currentHour*60*60) + (currentMin*60) + (currentSec);
        
        totalnumberOfSecs = numberOfSecsOfDayInEarlierTime + numberOfSecsOfDayInLaterTime;
        
    case 'present'
        totalnumberOfSecs = 0;
        
    case 'future'
        % currentYear:YearRef
        numberOfSecsOfDayInEarlierTime = totalNumberOfSecsIn1Day - (currentHour*60*60) - (currentMin*60) - currentSec;
        numberOfSecsOfDayInLaterTime = (HourRef*60*60) + (MinuteRef*60) + (SecondRef);
        
        totalnumberOfSecs = numberOfSecsOfDayInEarlierTime + numberOfSecsOfDayInLaterTime;    
end

end

%{
Note:
    The logic subsumes the time interval by 
        yearEnd, monthEnd, dayEnd - timeEarlier (within that year)
        + time interval (excluding that 2 years)
        + timeLater - yearEnd, monthEnd, dayEnd (within that year)

* 'This' refers earlier time
* 'That' refers later time

Caveat: does not include leap secs
%}
