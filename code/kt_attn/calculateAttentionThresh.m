clear;
close all;

load('subjects.mat');
N = length(subjectids);

% TODO: should calculate thresh per student?
all_attention = [];

for i=1:N
    sid=subjectids{i};
    filename = char(strcat('subjects/', sid, '.mat'));
    load(filename, 'data');
    
    haseeg = boolean(data.haseeg);
    attention=double(data.attention(haseeg));
    all_attention = [all_attention; attention];
    
    if (any(attention>90))
        sid
    end
end

p=mean(all_attention)
stddev=std(all_attention)

figure; hold on;
plot(all_attention, 'o');
plot([1 length(all_attention)], [p p]);

attentionThresh = p;
save('attentionThresh', 'attentionThresh');