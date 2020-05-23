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
%p = gcp();
%poolsize = p.NumWorkers;
poolsize = 8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%ESTABLISH PATHS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
addpath(fullfile(basepath,'src_fun', 'OPTIMIZATION'))
addpath(fullfile(basepath,'src_fun','STRUCTURE'))



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

%insert data
params.poolsize = poolsize;
params.compute_true_logZ = logical(params.compute_true_logZ);
params.data = data;
params.UDF = UDF;
params.coords = coords;
params.name = name;
params.data_directory = data_directory;
params.Filename = filename;
params.source_directory = source_directory;
params.exptdir=exptdir;

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
best_model = struct(best_model);

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

%I define this DeltaPEi using the elastic potential energy equation, where
%DeltaPE = kx^2 [ (vx)kx ] where k is a synaptic weight matrix and x is the state displacement. Thus, the frustration of any neuron in the
%system is described by the displacement of a spring, which will flow in
%the next timestep. 

%One can easily use a standard deltaEi from hopfield models, where 
%Delta Ei is equal to -(delta)Si*sum(synaptic weight matrix * Sj)

%It is not clear to me which is better yet. 

%DeltaPEi is the change in potential energy
p=0.2;
[PCNs] = PAPS(best_model,results,params,p);

%NOW WE ANALYZE THE GLOBAL ENERGY LANDSCAPE TO LABEL ATTRACTOR BASINS WITH
%THEIR ENSEMBLE COUNTERPARTS
load(params.Filename);
[SIM_STRUC] = simulation_analysis(best_model,results,params,STATE_TEMP,PCNs,p);

%[GE] = global_energy(best_model,data);

%VISUALIZE IT



%figure
%[C,h] = contour(GE)
%h.LineWidth=3;
%h.Fill='on';
%hold on
%X = [7 8 8 8 8 9 9 9 10 10];
%Y = [10 9 8 7 6 5 4 3 2 1];
%M=plot(X,Y)
%M.LineWidth=5;
%hold off

%We can also do somethings with this spring formulation,
%quantifying how much energy is shunted into the ensemble, provided we
%correct for a neural circuit being an open-thermodynamic system in which
%(practically speaking) energy is not constant.
save(strcat(params.exptdir, '/', 'results.mat'));