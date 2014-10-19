%
% Set the paths to several subdirectories within the core with
% respect to the main one.
%

function [] = setCorePaths()

information = localInformation('justTheFactsMam');
corePath    = information.codeplace;

% code examples
addpath(sprintf('%s/Examples',corePath));

% plotting code
addpath(sprintf('%s/Plots',corePath));

% useful functions shared by a lot of code, IDM transformations
addpath(sprintf('%s/Utils',corePath));

% neural network code
addpath(sprintf('%s/Netlab',corePath));

% SVM code
addpath(sprintf('%s/SVM',corePath));

% linear decompositions code - upcoming
%addpath(sprintf('%s/LinearDecompositions',corePath));
%addpath(sprintf('%s/LinearDecompositions/FastICA',corePath));

% old classifiers
addpath(sprintf('%s/Old_SVM',corePath));
addpath(sprintf('%s/Old_nnets',corePath));
