global DEBUG CACHE VERBOSE INDEX;
CACHE = expt.cache;
VERBOSE = '-printTicToc';
INDEX = {};

%% Setup
% get path of base folder
PATH = sprintf('%s/../', fileparts(mfilename('fullpath')));

% load libs
addpath(sprintf('%s/lib/cache', PATH));
addpath(sprintf('%s/lib/classifier', PATH));
addpath(sprintf('%s/lib/io', PATH));
addpath(sprintf('%s/lib/plot', PATH));
addpath(sprintf('%s/lib/stat', PATH));
addpath(sprintf('%s/lib/string', PATH));
addpath(sprintf('%s/lib/util', PATH));

addpath(sprintf('%s/ext/Calculate_Time_Interval', PATH));
addpath(sprintf('%s/ext/fdr_bh', PATH));
addpath(sprintf('%s/ext/fmri_core_new', PATH));
addpath(sprintf('%s/ext/hashtable', PATH));
addpath(sprintf('%s/ext/kappa', PATH));
addpath(sprintf('%s/ext/mult_comp_perm_t1', PATH));

warning('off', 'MATLAB:polyfit:PolyNotUnique');

seed = RandStream('mt19937ar','Seed', 5489);
RandStream.setGlobalStream(seed);

if exist('VERBOSE', 'var') && ~isempty(strfind(VERBOSE, '-printTicToc')),
	tic;
end

% load defaults
%expt.result = '../result';
if ~isfield(expt, 'desc'),
    expt.desc = sprintf('DESC %s', 'stuff'); %TODO figure out what to put in here
end;
expt

% make necessary folders
if ~exist(expt.cache, 'dir'),
    mkdir(expt.cache);
end;
if isfield(expt, 'classifier_folder') && ~exist(expt.classifier_folder, 'dir')
    mkdir(expt.classifier_folder);
end

% clear experimental cache from previous runs
clear_experimental_cache();
