function hasPassed = leapDayHasPassed( month, day)
%LEAPDAYHASPASSED assumes that the year is a leap year
% hasPassed = {0, 1} = {false, true}

if (month < 2)
    hasPassed = 0;
    return;
elseif (month > 2)
    hasPassed = 1;
    return;
else % (month) == feb
   hasPassed = 0;
   return;
%     if (day < 29)
%         hasPassed = 0;
%         return;
%     else % it is the leap day !
%         hasPassed = 0;
%         return;        
%     end
end

end

