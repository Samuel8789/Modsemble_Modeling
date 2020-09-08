%temp
%DONT RUN ME IN LINUX TERMINAL

warning('off','all');

basepath = pwd;

% THIRDPARTY DEPENDENCIES
addpath(fullfile(basepath,'thirdparty'))
addpath(fullfile(basepath,'thirdparty','QPBO-v1.32.src'))
addpath(fullfile(basepath,'thirdparty','glmnet_matlab'))
addpath(fullfile(basepath,'thirdparty','glmnet_matlab','glmnet_matlab'))
addpath(fullfile(basepath,'expt'))

% SOURCE FUNCTIONS

addpath(fullfile(basepath,'src_fun'))
addpath(fullfile(basepath,'src_fun','ANALYSIS'))
addpath(fullfile(basepath,'src_fun','FRAMEWORK'))
addpath(fullfile(basepath,'src_fun','MLE_STRUC'))
addpath(fullfile(basepath,'src_fun','include'))
addpath(fullfile(basepath,'src_fun','STRUCTURE'))

% d = uigetdir(pwd,'Select Data Folder');
d='/home/sp3660/Modsemble_Modeling/expt/darik_demo';
addpath(d);

load(strcat(d,'/','model_parameters.mat'));
load(strcat(d,'/','best_model.mat'));

%MUST HAVE PRELOADED EXT IDs

%Analysis for Core neurons//Pattern Completers

[ext_IDs] = find_IDs(params.data)
UDF = ext_IDs;

[~, results] = find_core_ext_IDs(best_model,params.data,UDF);

fprintf('Finding Pattern Completers')
fprintf('\n')

[PCNs, PAPS_INDEXED] = PAPS_score(best_model,results,params);
results.PCNs = PCNs;
results.PAPS_INDEXED = PAPS_INDEXED;
save(strcat(d,'/','results.mat'));

fprintf('\n')
fprintf('Core Analysis Completed')
fprintf('\n')

%temporary

temp_plot;