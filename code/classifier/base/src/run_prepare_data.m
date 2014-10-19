function [data] = run_prepare_data( expt )
%% Experiment

% load input files
task_data = read_task(expt.task_file);
eeg_data = read_eeg(expt.eeg_file, expt.sigqual_sample);

% prepare data
eeg_data = smooth_eeg(eeg_data, expt.denoise);

data = align_data(task_data, eeg_data);

data = gen_epochs(data, expt.epochs);
data = gen_epoch_features(data, expt.bands, expt.sampling_rate);
data = gen_higher_order_features(data, expt.higher_order_features);

data = gen_feature_matrix(data, expt.cond, expt.bands, expt.higher_order_features, expt.normalize);
data = filter_data(data, expt.sigqual_percentage);

end
