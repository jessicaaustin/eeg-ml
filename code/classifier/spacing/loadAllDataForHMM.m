

clear;
tic
[~, hmmData_2011_2012] = prepareData( '../../../data-lg/eeg_lex_view_2011_2012_cols.xls' );
save('hmmData_2011_2012');
toc; 

clear;
tic
[~, hmmData_2012_2013] = prepareData( '../../../data-lg/eeg_lex_view_2012_2013_cols.xls' );
save('hmmData_2012_2013');
toc; 

clear;
tic
[~, hmmData_2013_2014] = prepareData( '../../../data-lg/eeg_lex_view_2013_2014_cols.xls' );
save('hmmData_2013_2014');
toc

load gong.mat;soundsc(y);


% % group students across data sets
% load('hmmData_2011_2012');
% load('hmmData_2012_2013');
% load('hmmData_2013_2014');
% 
% hmmData = hmmData_2011_2012;
% for i=1:length(hmmData_2012_2013)
%     data = hmmData_2012_2013{i};
%     s = data.subjectID;
%     N = length(hmmData);
%     for j = 1:N
%         if strcmp(hmmData{j}.subjectID, s)
%             existingData = hmmData{j};
%             % TODO: should really group by word as well...
%             existingData.words = [existingData.words; data.words];
%             existingData.accept = [existingData.accept; data.accept];
%             existingData.latency = [existingData.latency; data.latency];
%             existingData.fluent = [existingData.fluent; data.fluent];
%             hmmData{j} = existingData;
%         else
%             hmmData{length(hmmData)+1} = data;
%         end
%     end
% end
