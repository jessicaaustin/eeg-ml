clear;
close all;

load('subjects.mat');
N = length(subjectids);

all_p = [];

for i=1:N
    sid=subjectids{i};
    fprintf('%s: ', char(sid));

    filename = char(strcat('subjects/', sid, '.mat'));
    load(filename, 'data');
    
    haseeg = boolean(data.haseeg);
    timeelapsed=data.timeelapsedms(haseeg);
    timeelapsed = timeelapsed/mean(double(timeelapsed));
    attention=data.attention(haseeg);
    attention = attention/mean(double(attention));
    meditation=data.meditation(haseeg);
    
%     figure; hold on;
%         plot(timeelapsed, attention, 'bo');
%         plot(timeelapsed, meditation, 'ro');
        p=double(timeelapsed)\double(attention);
%         plot(timeelapsed, p.*timeelapsed, 'g--');
%         ylim([min(attention), max(attention)]);
        fprintf('%0.2f\n', p);
        
        all_p = [all_p; p];
end

if false
    figure; hold on;
        plot(timeelapsed, attention, 'bo');
        p=double(timeelapsed)\double(attention);
        plot(timeelapsed, p.*timeelapsed, 'go--');
%         ylim([min(attention), max(attention)]);
        fprintf('%0.2f\n', p);
end

fprintf('Time elapsed and attention are correlated with r = %0.2f +/- %0.4f\n', ...
    mean(all_p), std(all_p));