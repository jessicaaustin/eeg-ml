function earlierTime = givesTheEarlierTime (time1, time2)
%TIMEDIFFERENCE 
currentMoment = time1;
T = time2;

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
    earlierTime = T;
elseif (elapsedTime == 0)   % isPresent
    earlierTime = currentMoment; % either will do
else
    earlierTime = currentMoment;
end


end

