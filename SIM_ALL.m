
%LOAD DIRECTORY
mod_dir = uigetdir(pwd,'SELECT_RESULT_DIRECTORY');
Files = dir(fullfile(mod_dir, '*.mat'));

MODEL_RESULTS = struct;
MODEL_RESULTS.SIMS = {};

for i = 1:length(Files)
    FILENAME = strcat(Files(i).folder, '/', Files(i).name);
    load(FILENAME);
    SIM_SCRIPT;
    MODEL_RESULTS.SIMS{end+1} = MODEL_STRUCTURE;
end
