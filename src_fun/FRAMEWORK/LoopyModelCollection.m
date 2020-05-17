classdef LoopyModelCollection

    properties
        %Properties & Definitions
        Num_Nodes, UDF_Count; %Neuronal + UDF Node Counts
        x_train, x_test; %Training and Test Data
        variable_names; %Node Names
        GLM_array; %Array of Learned Structure
        next_parameters; %Parameters for next Optimization Loop
        variable_groups; %Index of Possible Hyperedge Connections
        poolsize; %Parrallel Pool Size
        LASSO_options; %Structure Learning Settings
        s_lambda_sequence_LASSO; %Total Sequence of LASSO parameters
        s_lambda_sequence; %Iteration Specific Sequence of LASSO parameters
        p_lambda_sequence; %Iteration Specific Sequence of LASSO parameters
        density_sequence; %Archaic
        models;
        computed_true_logZ;%archaic
        hidden_model;%arhaic
        best_model_slambda, best_model_plambda; %Best Model Parameters
        L; %L
        next_s_parameters;
        next_p_parameters;
        learned_structures;
        logspace;
        value;
        fit_bois, gof_bois, graph_max;

    end

    methods
        function self = LoopyModelCollection(x_train, x_test, params, variable_groups)
            self.variable_groups = variable_groups;
            self.x_train = x_train;
            self.x_test = x_test;
            self.UDF_Count = size(params.UDF,2);
            self.Num_Nodes = size(params.data,2);
            self.variable_names = num2cell([1:size(x_train,2)]);
            self.density_sequence = [1];
            self.poolsize = params.poolsize;
            self.best_model_slambda = params.best_model_slambda;
            self.best_model_plambda = params.best_model_plambda;
            self.logspace=params.logspace;

            % lambdas are samples in a logspace
            if self.logspace == 1
                p_lambda_count = params.p_lambda_count;
                p_lambda_min_exp = log10(params.p_lambda_min);
                p_lambda_max_exp = log10(params.p_lambda_max);
                self.p_lambda_sequence = logspace(p_lambda_min_exp, p_lambda_max_exp, p_lambda_count);
                s_lambda_count = params.s_lambda_range;
                s_lambda_min_exp = log10(params.s_lambda_min);
                s_lambda_max_exp = log10(params.s_lambda_max);
                self.s_lambda_sequence_LASSO = logspace(s_lambda_min_exp, s_lambda_max_exp, s_lambda_count);
                self.s_lambda_sequence = params.s_lambda_min;
            else
                %lambdas are samples in a linspace
                p_lambda_count = params.p_lambda_count;
                s_lambda_count = params.s_lambda_range;
                self.p_lambda_sequence = linspace(params.p_lambda_min, params.p_lambda_max, p_lambda_count);
                self.s_lambda_sequence_LASSO = linspace(params.s_lambda_min, params.s_lambda_max, s_lambda_count);
            end
       
            %Build LASSO options
            opts.lambda = self.s_lambda_sequence_LASSO;
            self.LASSO_options = glmnetSet(opts);
            self.models = {};
        end
       
        function self = learn_priming_structures(self)

fprintf('\nLearning Structures using LASSOs\n');
%Here we learn a structure for every parameter in our selected sequence
[self.GLM_array] = lasso_node_by_node_group(self.x_train, self.variable_groups, self.LASSO_options);

%Here we ideally select the most sparse lambda. The smallest lambda is
%selected as the first lambda after 20% of the sequence has been learned to account for
%imprecision during the warm start. We also select the largest lambda. This
%ensures sufficient bracketing for the optimization procedures. 
%A stochastic set of structures are further selected to match the poolsize.

fprintf('\n Selecting Priming Structure\n');

max_s = length(self.s_lambda_sequence_LASSO);
min_s = fix(0.2*length(self.s_lambda_sequence_LASSO));
max_s_lambda = self.s_lambda_sequence_LASSO(max_s);
min_s_lambda = self.s_lambda_sequence_LASSO(min_s);
self.s_lambda_sequence = [max_s_lambda min_s_lambda];

 Pl = self.poolsize;
 selection = linspace(self.s_lambda_sequence_LASSO(min_s+1),self.s_lambda_sequence_LASSO(max_s-1),(Pl-2));
 self.s_lambda_sequence = [self.s_lambda_sequence selection];
 learned_structures = cell(1,self.poolsize);
 
  parfor i = 1:(Pl)
               learned_structures{i} = learn_structures(self.x_train, self.s_lambda_sequence(i), self.variable_groups, self.GLM_array,self.UDF_Count);
  end
  
  %PRE-ALLOCATING MODELS
              for j = 1:numel(self.s_lambda_sequence)
                  for k = 1:numel(self.p_lambda_sequence)
                      model = struct();
                      model.s_lambda = self.s_lambda_sequence(j);
                       model.p_lambda = self.p_lambda_sequence(k);
                       model.density = 1;
                        model.structure = learned_structures{j};
                        model.max_degree = max(sum(model.structure));
                        model.median_degree = median(sum(model.structure));
                        model.mean_degree = mean(sum(model.structure));
                        model.rms_degree = rms(sum(model.structure));
                       model.pending_parameter_estimation = true;
                       self.models{end+1} = model;
                   end
              end
                
              %REMOVE STRUCTURES THAT WERE TOO COLLINEAR TO INDUCE SPARSITY
              %SUCH THAT THE REGULARIZATION FAILS AND MAX DEGREE = 0
              for i = 1:length(self.models)
                  if self.models{i}.max_degree <= 0
                      self.models{:,i}={};
                  end
              end
              
             self.models = self.models(~cellfun('isempty',self.models));
              
              %REPLACE REMOVED STRUCTURES
              if length(self.models) < (Pl*2)
                  P=(((Pl*2)-(length(self.models)))/2);
                  selection = randsample(self.s_lambda_sequence_LASSO((min_s+1):(max_s-1)),P);
                  L = numel(self.s_lambda_sequence);
                  self.s_lambda_sequence = [self.s_lambda_sequence selection];
                  for i = 1:P
                      learned_structures{L+i} = learn_structures(self.x_train, selection(i), self.variable_groups, self.GLM_array,self.UDF_Count);
                  end
                  
                  %REPLACE REMOVED MODELS
                  
                  for j = (L+1):(L+P)
                      for k = 1:numel(self.p_lambda_sequence)
                          model = struct();
                          model.s_lambda = self.s_lambda_sequence(j);
                          model.p_lambda = self.p_lambda_sequence(k);
                          model.density = 1;
                          model.structure = learned_structures{j};
                          model.max_degree = max(sum(model.structure));
                          model.median_degree = median(sum(model.structure));
                          model.mean_degree = mean(sum(model.structure));
                          model.rms_degree = rms(sum(model.structure));
                          model.pending_parameter_estimation = true;
                          self.models{end+1} = model;
                   end
                  end
              else
                  L = numel(self.s_lambda_sequence);
                  P = 0;
              end
              self.L = L+P;
              self.learned_structures=learned_structures;
        end
        
        function self = do_loopy_structure_learning_shuffle(self)
            fprintf('Structure learning using Lasso Logistic Regression\n');
            
            
            %specific parameters for shuffled models
            
            fprintf('\nSelecting Priming Structures\n')
            
            self.learned_structures = ones(size(self.x_train,2),size(self.x_train,2));
           
                        model = struct();
                        model.s_lambda = self.best_model_slambda;
                        model.p_lambda = self.best_model_plambda;
                        model.density = 1;
                        model.structure = self.learned_structures;
                        model.max_degree = max(sum(model.structure));
                        model.median_degree = median(sum(model.structure));
                        model.mean_degree = mean(sum(model.structure));
                        model.rms_degree = rms(sum(model.structure));
                        model.pending_parameter_estimation = true;
                        self.models{end+1} = model;
                    
                end
                
        function self = do_parameter_estimation(self, BCFW_max_iterations, BCFW_fval_epsilon,...
                computed_true_logZ, reweight_denominator, printInterval, printTest, MaxTime)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%PARAMETER ESTIMATION%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

           %Iterate through non-modeled models and estimate parameters
           %(lol, tongue twister)
            for i = 1:numel(self.models)
                if self.models{i}.pending_parameter_estimation == true
                    model = self.models{i};
                    
                    fprintf('\nParameter Estimation: s_lambda=%e; p_lambda=%e\n',...
                    model.s_lambda, model.p_lambda);
               
                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%OVERPARAMETERIZATION%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               

                % First we need to overcomplete the structure/samples to
                % run parameter estimation. 
                
                % We can save time by using an overcomplete_struct from the
                % prior model if the s_lambda is identical
                
                
               
                if i == 1 || any(any(self.models{i-1}.structure ~= model.structure))
                    overcomplete_struct = samples_to_overcomplete(self.x_train, model.structure);
                end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%REWEIGHTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Next we need to reweight to ensure a convex free-energy

                % define the reweight parameter
                if ischar(reweight_denominator)
                    if strcmp(reweight_denominator, 'max_degree')
                        reweight = 2/model.max_degree;
                    elseif strcmp(reweight_denominator, 'median_degree')
                        reweight = 2/model.median_degree;
                    elseif strcmp(reweight_denominator, 'mean_degree')
                        reweight = 2/model.mean_degree;
                    elseif strcmp(reweight_denominator, 'rms_degree')
                        reweight = 2/model.rms_degree;
                    else
                        error('Unknown reweighting denominator ''%s''', reweight_denominator);
                    end
                else
                   reweight = 2/reweight_denominator;
                end

                % sanity check
                if reweight > 1
                    reweight = 1;
                end

                % saves the reweight for future cross-validation
                model.reweight = reweight;

               % Now we create our model object.
                loopy_model_train_object = Ising( ...
                    overcomplete_struct.YN, ...
                    overcomplete_struct.YE, ...
                    overcomplete_struct.Ut, ...
                    overcomplete_struct.Vt, ...
                    overcomplete_struct.Ns, ...
                    overcomplete_struct.edges, ...
                    model.p_lambda, ...
                    'checkStuck', false, ...
                    'reweight', reweight);
             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%BCFW%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Now we actually do the meat n potatoes

%First we create/run the BCFW on our object

                bcfw = BCFW(loopy_model_train_object, ...
                    'printInterval', printInterval, ...
                    'printTest', printTest, ...
                    'printComputedDualityGap',true,...
                    'MaxTime',MaxTime,...
                    'MaxIter', BCFW_max_iterations, ...
                    'fvalEpsilon', BCFW_fval_epsilon);
                
                bcfw.run();

%We now extract F and G, which score the two neuron states (ON/OFF) and
%four connectivity states (ON-OFF, ON-OFF, OFF-ON, OFF-OFF). These are in
%the format of 2xNode 4xEdge by the way. These are the Phi's. 

                model.theta = bcfw.obj.computeParams();

%We now approx. the partion function which globally
%normalizes our scores
                
                logZ = bcfw.obj.partition_function(model.theta);

%We now convert F and G into the node and edge potentials of our model, and
%adjust partition function. Now we have a nice graphical model.
                fprintf('Converting F and G to node and edge potentials\n');
                [node_pot, edge_pot, logZ_pot] = get_node_and_edge_potentials(model.theta.F,...
                    model.theta.G, logZ, overcomplete_struct.edges{1}');
                model.theta.node_potentials = node_pot;
                model.theta.edge_potentials = edge_pot;
                model.theta.logZ = logZ_pot;

%If you want you can compute exact true_logZ and true_node_marginals, but
%this is very computationally expensive and you almost certainly shouldn't.

                if self.computed_true_logZ
                    fprintf('Starting to run JTA to compute true partition function\n');
                    [true_node_marginals,~,~,~,~,~,true_logZ] = run_junction_tree(node_pot, edge_pot, 'verbose', true);
                    model.theta.true_logZ = true_logZ;
                    model.true_node_marginals = true_node_marginals;
                end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%ASSESS MODEL%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
                
                    % Compute training likelihood
                    fprintf('Computing training likelihood\n');
                    model.train_likelihood = compute_avg_log_likelihood( ...
                        model.theta.node_potentials, ...
                        model.theta.edge_potentials, ...
                        model.theta.logZ, ...
                        self.x_train);

                    % Compute test likelihood
                    fprintf('Computing test likelihood\n');
                    model.test_likelihood = compute_avg_log_likelihood( ...
                        model.theta.node_potentials, ...
                        model.theta.edge_potentials, ...
                        model.theta.logZ, ...
                        self.x_test);
               
                    model.pending_parameter_estimation = false;
                    self.models{i} = model;
                end
            end
            fprintf('Finished estimating parameters.\n');
        end
        
        function self = learn_optimized_structures(self)
            fprintf('\nLearning Structures Using Optimized Parameters\n');
            P = numel(self.next_s_parameters);
            PP = numel(self.next_p_parameters);
            learned_structures = self.learned_structures;
            for i = self.L+1:(self.L+P)
                learned_structures{i} = learn_structures(self.x_train, self.s_lambda_sequence(i) ,self.variable_groups, self.GLM_array,self.UDF_Count);
            end
            
             % we will initialize structures for all new. This means we will
             % replicate the newly learned structures |old (step_1) and new(step_2) p_lambda_sequence| times
             
             %STEP1
            if ismember(self.next_p_parameters,self.p_lambda_sequence) ~= 1
                for j = 1:self.L
                    for k = 1:PP
                        model = struct();
                        model.s_lambda = self.s_lambda_sequence(j);
                        model.p_lambda = self.next_p_parameters(k);
                        model.density = 1;
                        model.structure = learned_structures{j};
                        model.max_degree = max(sum(model.structure));
                        model.median_degree = median(sum(model.structure));
                        model.mean_degree = mean(sum(model.structure));
                        model.rms_degree = rms(sum(model.structure));
                        model.pending_parameter_estimation = true;
                        self.models{end+1} = model;
                    end
                end
                self.p_lambda_sequence = [self.p_lambda_sequence self.next_p_parameters];
            end
            
               %STEP2
                for j = self.L+1:self.L+P
                    for k = 1:length(self.p_lambda_sequence)
                        model = struct();
                        model.s_lambda = self.s_lambda_sequence(j);
                        model.p_lambda = self.p_lambda_sequence(k);
                        model.density = 1;
                        model.structure = learned_structures{j};
                        model.max_degree = max(sum(model.structure));
                        model.median_degree = median(sum(model.structure));
                        model.mean_degree = mean(sum(model.structure));
                        model.rms_degree = rms(sum(model.structure));
                        model.pending_parameter_estimation = true;
                        self.models{end+1} = model;
                    end
                end
                
                 %REMOVE STRUCTURES THAT WERE TOO COLLINEAR TO INDUCE SPARSITY
              %SUCH THAT THE REGULARIZATION FAILS AND MAX DEGREE = 0
              for i = 1:length(self.models)
                  if self.models{i}.max_degree <= 0
                      self.models{:,i}={};
                  end
              end
              
             self.models = self.models(~cellfun('isempty',self.models));
             self.L = self.L+P;
             self.learned_structures=learned_structures;
        end

        function self =  maximimum_likelihood_optimize(self)
            L = numel(self.models);
            pb =[];
            sb=[];
            ttb=[];
            trnb=[];
            for i = 1:L
                pb = [pb self.models{i}.p_lambda];
                sb = [sb self.models{i}.s_lambda];
                ttb = [ttb self.models{i}.test_likelihood];
                trnb = [trnb self.models{i}.train_likelihood];
                i=i+1;
            end
            %difference in prediction
            LL_diff = abs(trnb-ttb);
            %min-max normalized (LL_diff max-min normalized)
            LL = (trnb-min(trnb))/(max(trnb)-min(trnb)).*(3/6)+(ttb-min(ttb))/(max(ttb)-min(ttb)).*(2/6)+(LL_diff-max(LL_diff))/(min(LL_diff)-max(LL_diff)).*(1/6);
         
            
%LL = transpose(LL);
Z = [sb;pb;LL];
Z = transpose(Z);
unique_p_cell = numel(unique(pb));
unique_s_cell = numel(unique(sb));
Z = sortrows(Z,[1 2]);
X = linspace(unique_p_cell,unique_p_cell,unique_s_cell);
ZZ = mat2cell(Z,X);
L3 = numel(ZZ);
fit_bois = cell(1,L3);
gof_bois = cell(1,L3);
for m = 1:L3
    p = ZZ{m}(:,2);
    l = ZZ{m}(:,3);
    [fitresult, gof] = createFit_C(p,l,0)
    fit_bois{m} = fitresult;
    gof_bois{m} = gof;
end

Z2 = sortrows(Z,[2 1]);
Y = linspace(unique_s_cell,unique_s_cell,unique_p_cell);
ZZZ = mat2cell(Z,Y)

L4 = numel(ZZZ)
for m = 1:L4
    s = ZZZ{m}(:,1);
    l = ZZZ{m}(:,3);
    [fitresult, gof] = createFit_B(s,l,0)
    fit_bois{end+1} = fitresult;
    gof_bois{end+1} = gof;
end
p_max = [];
s_max = [];
g_max = [];
graph_max = [];
P_X = linspace(1,max(pb),1000);
S_X = linspace(1,max(sb),1000);
usb = unique(sb);
upb = unique(pb);
for n = 1:numel(fit_bois)
    if n <=L3
        [max_idx,g_max] = find_graph_cut_max(fit_bois{n},P_X);
        p_max = [p_max max_idx];
        s_max = [s_max usb(n)];
        graph_max = [graph_max g_max];
    else n > L3
        [max_idx,g_max] = find_graph_cut_max(fit_bois{n},S_X);
        s_max = [s_max max_idx];
        p_max = [p_max upb(n-numel(usb))];
        graph_max = [graph_max g_max];
    end
end

graph_max = [graph_max;s_max;p_max];
self.fit_bois=fit_bois;
self.gof_bois=gof_bois;
self.graph_max = graph_max;
        end
        
        function self = select_next_parameters(self)
            [M,I] = max(self.graph_max(1,:),[],2);
            self.value = M;
            self.next_s_parameters=self.graph_max(2,I);
            self.next_p_parameters=self.graph_max(3,I);
            self.s_lambda_sequence=[self.s_lambda_sequence self.next_s_parameters];
            
            %Now we add parameters from gaussian distributions centered
            %with a mean I [s,p] and  sigma equal 0.5(range(s,p))/7;
            
            self.next_s_parameters = [self.next_s_parameters normrnd(self.graph_max(2,I),(0.5*(params.s_lambda_max-params.s_lambda_min))/7,[1 floor(sqrt(params.poolsize*2))])];
            self.next_p_parameters = [self.next_p_parameters normrnd(self.graph_max(3,I),(0.5*(params.p_lambda_max-params.p_lambda_min))/7,[1 floor(sqrt(params.poolsize*2))])];
            
        end
        
    end
end

