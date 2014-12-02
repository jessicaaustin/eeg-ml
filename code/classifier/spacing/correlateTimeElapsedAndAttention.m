clear;
close all;

load('subjects.mat');
N = length(subjectids);


for i=1:N
    sid=subjectids{i};
    fprintf('%s: ', char(sid));

    filename = char(strcat('subjects/', sid, '.mat'));
    load(filename, 'data');
    
    haseeg = boolean(data.haseeg);
    timeelapsed=data.timeelapsedms(haseeg);
    attention=data.attention(haseeg);
    meditation=data.meditation(haseeg);
    
%     figure; hold on;
%         plot(timeelapsed, attention, 'bo');
%         plot(timeelapsed, meditation, 'ro');
        p=double(timeelapsed)\double(attention);
%         plot(timeelapsed, p.*timeelapsed, 'g--');
%         ylim([min(attention), max(attention)]);
        fprintf('%0.2f\n', p);
end