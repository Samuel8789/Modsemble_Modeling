%% MODSEMBLE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODSEMBLE (temp name until sassy acronym)
% See README for papers and proper citation
% Laboratory of Rafael Yuste
% Neurotechnology Center
% Columbia University, NY, NY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% STARTUP

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%CREATE PARALLEL POOLS%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create parallel pool
p = gcp();
poolsize = p.NumWorkers;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%ESTABLISH PATHS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

basepath = pwd;

% THIRDPARTY DEPENDENCIES
addpath(fullfile(basepath,'thirdparty'))
addpath(fullfile(basepath,'thirdparty','QPBO-v1.32.src'))
addpath(fullfile(basepath,'thirdparty','glmnet_matlab'))
addpath(fullfile(basepath,'thirdparty','glmnet_matlab','glmnet_matlab'))

% SOURCE FUNCTIONS

addpath(fullfile(basepath,'src'))
addpath(fullfile(basepath,'src','ANALYSIS'))
addpath(fullfile(basepath,'src','FRAMEWORK'))
addpath(fullfile(basepath,'src','MLE_STRUC'))
addpath(fullfile(basepath,'src','include'))
addpath(fullfile(basepath,'src', 'OPTIMIZATION'))
addpath(fullfile(basepath,'src','STRUCTURE'))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%DATA IMPORT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Select File
filename = uigetfile('*.mat','Select Data File');
load(filename);
[pathstr,name,ext] = fileparts(filename);

%Data Directory
data_directory = uigetdir(pwd,'Select Data Directory');
data_directory = strcat(data_directory, '/');

%Source Directory
source_directory = uigetdir(pwd,'Select Modsemble Directory');
source_directory = strcat(source_directory, '/');

%Create Results Folders
exptdir = strcat(source_directory, 'expt', '/', name);
temp_folder = strcat(source_directory, 'expt', '/', name, '/tmp');
mkdir(temp_folder);
addpath(exptdir);
addpath(temp_folder);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%INITIALIZE PARAMETERS%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%MUST CHANGE MANUALLY
%load parameter structure
load('parameters.mat');

%% CONDITIONAL RANDOM FIELD MODELING

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%CRF MODELING%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Models] = run(params,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%SHUFFLE MODELS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%LOAD PARAMETERS OF BEST MODELS
load(strcat(params.exptdir,'/','best_model.mat'));
params.best_model_slambda = best_model.s_lambda;
params.best_model_plambda = best_model.p_lambda;

%RUN SHUFFLES IN PARALLEL
parfor i=1:params.number_of_shuffles
   run(params,1);
end

%IMPORT MODELS FROM TEMP DIRECTORY
shuffled_models = {};
shuffle_directory = uigetdir(pwd,'SELECT_TEMP_DIRECTORY');
Files = dir(fullfile(shuffle_directory, '*.mat'));
for i = 1:length(Files)
    FILENAME = strcat(Files(i).folder, '/', Files(i).name);
    shuffled_MODEL = load(FILENAME);
    shuffled_models{end+1} = shuffled_MODEL.model_collection.models{1,1};
end

%SAVE AS ONE COLLECTION
save_shuffled_models(shuffled_models,params);
shuffle_model = load(strcat(params.source_directory,'/expt/', params.name, '/', 'fulldata.mat'));

%% FIND PATTERN COMPLETORS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%ATTRACTOR ANALYSIS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%ENSEMBLE ANALYSIS
clearvars -except params best_model shuffle_model
[results] = find_plot_temporal_crf_core(best_model,shuffle_model,params.data,params.UDF, params.coords);
save(strcat(params.exptdir, '/', 'results.mat'));

%The former gives some insight into pattern-completing neurons
%Now we do more dedicated analysis into them

%PATTERN COMPLETOR ANALYSIS
%I FIND TARGETS USING THE 0-1 PAPS SCALE, WHERE PAPS IS DEFINED AS
% PAPS = (n_NodeStr+n_Degree+n_deltaPEi+AUC)/4, where n_ designated min-max
% normalization

%DeltaPEi is the change in potential energy

[PCNs] = PAPS(best_model,results,params);

%NOW WE ANALYZE THE GLOBAL ENERGY LANDSCAPE TO LABEL ATTRACTOR BASINS WITH
%THEIR ENSEMBLE COUNTERPARTS


[GE] = global_energy(best_model,data);

%VISUALIZE IT



figure
[C,h] = contour(GE)
h.LineWidth=3;
h.Fill='on';
hold on
X = [7 8 8 8 8 9 9 9 10 10];
Y = [10 9 8 7 6 5 4 3 2 1];
M=plot(X,Y)
M.LineWidth=5;
hold off

