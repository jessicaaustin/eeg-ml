Calculate Time Interval
--------------------------------------
main executing reference usage: usage_timeDifference.m 

Sample output:
The time interval between earierTime [09-Sep-1900 09:09:18] and laterTime [31-Oct-2010 23:33:08] is 
110 years, 79 days, 14 hours, 23 mins, 50.17 secs.


The objective is to compute the time interval taking into account the leap days that are subsumed within the time interval. 

Note:
    The logic subsumes the time interval by 
        yearEnd, monthEnd, dayEnd - timeEarlier (within that year)
        + time interval (excluding that 2 years)
        + timeLater - yearEnd, monthEnd, dayEnd (within that year)

It differentiates the earlier time and later time, ie. timeDifference (earlierTime, laterTime) gives the same result as timeDifference (laterTime, earlierTime).

leapDayHasPassed.m illustrates a possible approach to determine if a certain day of concern is already over at a given time.

Caveat: does not include leap secs 


If the reference demo has a more elegant presentation, please do not hesitate to suggest and send feedback to author.
Email: promethevx@yahoo.com.

Thank you.

Regards,
Michael Chan JT

-------------------------------------- EOF --------------------------------------