 close all
 
%%

clear;
load('timeElapsed');
N = length(timeElapsed);

bins = 0:1:11;
timeElapsedPerSubject=zeros([length(bins)-1 N]);
for i=1:N
   s = timeElapsed{i};
   s(s>=10) = 10;
   for j=1:length(bins)-1
       timeElapsedPerSubject(j,i) = sum(s>=bins(j) & s<bins(j+1));
   end
end

figure; bar(timeElapsedPerSubject);
set(gca,'YScale','log');
ylabel('Count for each time (log scale)');
xlabel('Time elapsed between words (seconds)')


%% 

clear;
load('testSubject_vs_words');

figure;
[n, xout] = hist(subjectVSwords',[2:10]);
bar(xout, n, 'barwidth', 1, 'basevalue', 1);
set(gca,'YScale','log');
ylabel('Number of sequences (log scale)');
xlabel('Length of word sequence')