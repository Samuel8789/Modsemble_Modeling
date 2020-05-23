function [SIM_STRUC] = simulation_analysis(best_model,results,params,STATE_TEMP,PCNs,p)

% FIRST INITIALIZE PARAMETERS AND COUNTS
data=params.data;
UDF=params.UDF;
num_neur = size(data,2);
num_ens = size(UDF,2);
ens = results.core_crf;
ENS_STATE = cell(1,num_ens);
%STATE_TEMP=STATE_TEMP;
num_neur = size(data,2);
num_ens = size(UDF,2);
ens = results.core_crf;
ENS_STATE = cell(1,num_ens);
auc=results.auc(1:num_neur,:);

name = strcat(params.name,'.hdf5');

weights = h5read(name, '/weights');
fweights = weights.*best_model.structure(1:num_neur,1:num_neur);
pos_fweights = fweights;
edge_potentials = fweights;
pos_fweights(pos_fweights < 0)=0;
node_str = sum(pos_fweights);
degrees = sum(best_model.structure(1:num_neur,1:num_neur));

%STATE_TEMP = h5read('MODEL_LARGE.hdf5','/STATE_TEMP');

%find ON neurons in ensemble
for i = 1:num_ens
    ENS_STATE{i} = auc(:,i);
    ENS_STATE{i} = round(ENS_STATE{i});
    %ENS_STATE{i}=STATE_TEMP(:,6,i);
    %ENS_STATE{i}(ENS_STATE{i}==-1)=0;
    %ENS_STATE{i}=ENS_STATE{i}.*transpose([1:400]);
    %ENS_STATE{i}(ENS_STATE{i}==0)=[];
end

size_ens = cell(1,num_ens);
for i = 1:num_ens
    size_ens{i}=sum(round(auc(:,i)));
end

%transition times T
T=cell(1,31);
for i = 1:num_ens
    for ii = 2:size(data,1)
        if UDF(ii,i)==1 & UDF(ii-1,i)==0
            T{i}=[T{i} ii];
        end
    end
end

data=transpose(data);
data(data<1)=-1;

delta_Ei_normal = @(Si,t) -(data(Si,t)-data(Si,t-1))*sum(transpose(edge_potentials(Si,:)).*data(Si,t));

DELTA_U=cell(1,num_ens);

for i = 1:num_ens
    DELTA_U{i}=zeros(length(T{i}),num_neur);
    z=1;
    for ii = T{i}
        for k = 1:num_neur
            DELTA_U{i}(z,k) = delta_Ei_normal(k,ii);
        end
        z=z+1;
    end
end

for i = 1:num_ens
    ENS_STATE{i}=ENS_STATE{i}.*transpose([1:num_neur]);
    ENS_STATE{i}(ENS_STATE{i}==0)=[];
end

DELTA_ENS_SUM = cell(1,num_ens);
DELTA_U_SUM=cell(1,num_ens);
for i = 1:num_ens
    %for ii = 1:num_neur
    DELTA_U_SUM{i} = sum(DELTA_U{i});
    DELTA_ENS{i} = DELTA_U{i}(:,[ENS_STATE{i}]);
    DELTA_ENS_SUM{i}=sum(DELTA_ENS{i});
end

DEMA = []; %Max Ei
DEMI = []; %Min Ei
nE = []; %min-max norm Ei
dmax = []; %max deg
dmin = []; %min deg
nD = []; %min-max norm deg
nsmi = []; %min node str
nsma = []; %max node str
nS = []; %min max norm node str
auc2 = [];%auc by ens
PAPS_INDEXED = cell(2,num_ens);

for i = 1:num_ens
    DEMA = min(DELTA_ENS_SUM{i})*-1;
    DEMI = max(DELTA_ENS_SUM{i})*-1;
    nE = ((DELTA_U_SUM{i}([ENS_STATE{i}]).*-1)-DEMI)./(DEMA-DEMI);
    
    dmax = max(degrees([ENS_STATE{i}]));
    dmin = min(degrees([ENS_STATE{i}]));
    nD = (degrees(ENS_STATE{i})-dmin)./(dmax-dmin);
    
    nsmi= min(node_str([ENS_STATE{i}]));
    nsma = max(node_str([ENS_STATE{i}]));
    nS = transpose((node_str([ENS_STATE{i}])-nsmi)./(nsma-nsmi));
        
    auc2 = auc([ENS_STATE{i}],i);
    auc2=transpose(auc2);
    
    PAPS_INDEXED{1,i}=ENS_STATE{i};
    for ii = 1:length(ENS_STATE{i})
        PAPS = @(Ni) (nE(Ni)+nD(Ni)+nS(Ni)+auc2(Ni))/4;
        PAPS_INDEXED{2,i}=[PAPS_INDEXED{2,i} PAPS(ii)];
    end
end
   
SCNs = cell(1,num_ens);
for i = 1:num_ens
    [V,idx] = maxk(PAPS_INDEXED{2,i},ceil(p*size_ens{i}));
    SCNs{i} = transpose(PAPS_INDEXED{1,i}(idx));
end

SELECTIONS = {};

for i = 1:length(SCNs)
    SELECTIONS{end+1} = [SCNs{i}; PCNs{i}]; 
end

JIP = {};
for i = 1:length(SELECTIONS)
    JIP{end+1} = pdist((SELECTIONS{i}),'jaccard');
end

JI = [];
for i = 1:length(JIP)
    JI = [JI JIP{i}];
end
 
%JI = cell2mat(JI);
JIm = mean(JI);
JIstd = std(JI);
JIsem = JIstd/sqrt(length(JI));

SIM_STRUC.JIP = JIP;
SIM_STRUC.JI = JI;
SIM_STRUC.JIm=JIm;
SIM_STRUC.JIsem=JIsem;
SIM_STRUC.JIstd=JIstd;
SIM_STRUC.SELECTIONS=SELECTIONS;
end
    

        


