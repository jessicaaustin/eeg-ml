clear;

addpath('../base/src')
global EXPERIMENTAL;
EXPERIMENTAL = {};

% data loading
% expt.cache = '';
expt.cache = '../../../../cache'; % the location of the cache folder
expt.task_file = '../../../data/task_lex_view.xls'; % the location of the task file

% create cache, etc
run_setup;

% prepare data
[data, words, sequences] = prepareData( expt );
