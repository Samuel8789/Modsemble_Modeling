%CHOOSE FILE
%filename = uigetfile();
%load(filename);
load(params.Filename);

MODEL_STRUCTURE = struct;

%These are % of total ens_size to test
values = [0.40 0.30 0.25 0.20 0.15 0.10 0.05 0.025 0.01 0.005];

%This is my equations DeltaEi = -DeltaSi times Weight*Si
MODEL_STRUCTURE.ME={};
for p = 1:length(values)
    [PCNs] = PAPS(best_model,results,params,values(p));
    [SIM_STRUC] = simulation_analysis(best_model,results,params,STATE_TEMP,PCNs,values(p));
    MODEL_STRUCTURE.ME{end+1} = SIM_STRUC;
end

%This is hopfield DeltaEi = -DeltaSi times Weight times Sj
MODEL_STRUCTURE.HOP={};
for p = 1:length(values)
    [PCNs] = PAPS_hop(best_model,results,params,values(p));
    [SIM_STRUC] = simulation_analysis_hop(best_model,results,params,STATE_TEMP,PCNs,values(p));
    MODEL_STRUCTURE.HOP{end+1} = SIM_STRUC;
end

%This is spring -(DeltaSi^2)times weight
MODEL_STRUCTURE.SPRING={};
for p = 1:length(values)
    [PCNs] = PAPS_spring(best_model,results,params,values(p));
    [SIM_STRUC] = simulation_analysis_spring(best_model,results,params,STATE_TEMP,PCNs,values(p));
    MODEL_STRUCTURE.SPRING{end+1} = SIM_STRUC;
end

%jaccard sim for each run
JI_ME = [];
JI_HOP = [];
JI_SPRING = [];
for i = 1:length(MODEL_STRUCTURE.ME)
    JI_ME = [JI_ME MODEL_STRUCTURE.ME{i}.JIm];
    JI_HOP = [JI_HOP MODEL_STRUCTURE.HOP{i}.JIm];
    JI_SPRING = [JI_SPRING MODEL_STRUCTURE.SPRING{i}.JIm];
end

%values = values.*200;
%figure
%A = plot(values,JI_ME);
%hold on
%B = plot(values,JI_HOP);
%C = plot(values,JI_SPRING);
%hold off

