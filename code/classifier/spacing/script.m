clear
addpath('../base/src')
global EXPERIMENTAL;
EXPERIMENTAL = {};

% data loading
% expt.cache = '';
expt.cache = '../../../../cache'; % the location of the cache folder
expt.task_file = '../../../data/task_lex_view.xls'; % the location of the task file
expt.eeg_file = '../../../data/eeg.xls'; % the location of the eeg file

% data info
expt.sigqual_sample = 0; % sigqual is from 0 (best) to 200 (worst).  Any eeg segment (row in eeg file) with sigqual higher than sigqual_sample will be filtered out
expt.sigqual_percentage = 1.0; % Any task segment (row in task file) with a percentage of filtered out eeg segments higher than sigqual_percentage will be filtered out
expt.denoise = 1; % denoise by wavelet transform
expt.bands = [1 4 8 12 30 101]; % the boundaries of the frequency band buckets.  e.g. 1-4 is a bucket, 4-8 is the next bucket
expt.sampling_rate = 512; % sampling rate of the eeg device used to collect data
expt.normalize = 1; % per-subject normalization (by zscore) across entire dataset prior to cross-validation and balancing
expt.cond = 'cond'; % the name of the task file column that represents the dependent variable of the experiment

% segmental features
% each data segment is broken into sub-segments (epochs) and higher-order features are computed over these epochs 
expt.epochs.length = 1.0; % the length of an epoch (in seconds)
expt.epochs.overlap = 0.5; % the overlap between two epochs (in seconds)
expt.epochs.min = 0.5; % the minimum length of an epoch (in seconds)
expt.higher_order_features = {'mean', 'var'}; % mean, var (variance), max, and min are supported

% feature selection
expt.feature_selector.algorithm = 'pca'; % The feature selection algorithm.  Possible values: rank or pca
expt.feature_selector.dimensions = 3; % The dimensions to use after performing pca.  Field required if feature_selector.algorithm is pca
%expt.feature_selector.num_selected = 10; % The number of selected features in rank feature selection.  Field required if feature_selector.algorithm is rank

% classify
expt.balance = 'undersample'; % how to balanace the data.  Possible values: undersample, oversample, none
expt.cv = 'leave-one-out'; % the type of cross-validation.  Leave-one-out is the only available method currently
expt.cv_subjects = 'within'; % the type of train/test split: within subject or between subject
expt.classifier = 'svm'; % the type of algorithm.  Possible values: svm, nbayesPooled

% run the experiment
[data, results] = run_experiment(expt);
