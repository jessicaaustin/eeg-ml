function [ upperUnit, lowerUnit ] = adjustUnits2(upperUnit, lowerUnit, base)
%ADJUSTUNITS2 

while (lowerUnit >= base) 
    lowerUnit = lowerUnit - base;
    upperUnit = upperUnit + 1;
end

end

