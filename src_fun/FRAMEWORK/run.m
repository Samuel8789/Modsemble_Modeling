function [model_collection] = run(params,shuffle_status)
%% INITIALIZE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%INITIALIZE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% IMPORT DATA
        %X = Samples x Neuron Matrix
        X = params.data;
        
        %Length of Y
        Num_Samples = size(X,1);
        
        %Number of Neuronal Nodes
        Num_Nodes = size(X,2);
        
        %SEPARATE INTO TRAINING AND TEST SETS (WITHOLD FOR VALIDATION)
        x_train = X(1:floor(params.split*Num_Samples),:);
        x_test = X((floor(params.split*Num_Samples)+1):Num_Samples,:);
        
        %Determine whether there exists user-defined features (UDF)
        UDF_Count = size(params.UDF, 2);
        
        %Merge UDF and Neuronal Nodes
        if UDF_Count > 0
            assert(Num_Samples == size(params.UDF, 1), ...
                   'UDF and neuron data must have same number of samples.')
            x_train = [x_train params.UDF(1:floor(params.split*Num_Samples),:)];
            x_test = [x_test params.UDF((floor(params.split*Num_Samples)+1):Num_Samples,:)];
        end
        
 % SHUFFLE IF NECESSARY
 
        %Simple shuffle if generating shuffled models
        if shuffle_status == 1
            for i = 1:size(x_train,2);
                x_train(:,i) = randsample(x_train(:,i),size(x_train,1));
            end
        end
        
 % CONSTRAIN POTENTIAL EDGES
    
        %ALL EDGES ARE POSSIBLE
        variable_groups = all_but_me(1, size(x_train, 2));
        
        %UDF CANNOT CONNECT TO EACH OTHER, ONLY NEURONS
        if params.hyperedge == 2;
            for ii = ((Num_Nodes+1):(Num_Nodes+UDF_Count))
                variable_groups{ii} = [1:Num_Nodes];
            end
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%GENERATE MODELING OBJECT%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model_collection = LoopyModelCollection(x_train, x_test, params, variable_groups);        
      
%% STRUCTURE LEARNING
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%FIRST ROUND OF LEARNING%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if shuffle_status == 0     
    model_collection = model_collection.learn_priming_structures();
else shuffle_status == 1
    model_collection = model_collection.do_loopy_structure_learning_shuffle();
end

%% PARAMETER ESTIMATION

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%FIRST ROUND OF PARAEMETER ESTIMATION%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Parameter Estimation...\n');
%Training
        model_collection = model_collection.do_parameter_estimation(...
            params.BCFW_max_iterations, params.BCFW_fval_epsilon,...
            params.compute_true_logZ, params.reweight_denominator,...
            params.printInterval, params.printTest, params.MaxTime)
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%SAVE_RESULTS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
      if shuffle_status == 0
       save(strcat(params.exptdir, '/', 'model_collection.mat'), 'model_collection');
       [best_model_index] = get_best_model(model_collection);
       [best_model] = SingleLoopyModel(model_collection, best_model_index);
       save(strcat(params.exptdir, '/', 'best_model.mat'), 'best_model');
      else shuffle_status == 1
          save(strcat(params.exptdir, [tempname, '.mat']),'model_collection');
      end
      
%% OPTIMIZATION LOOPING

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%OPTIMIZE NEXT ITERATION%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if shuffle_status == 0
    if params.manual_iter > 0
        for i = 1:params.manual_iter
           % MANUAL ITER SET TO ENSURE REASONABLE RUNTIME
           
        model_collection =  model_collection.maximimum_likelihood_optimize()
        model_collection = model_collection.select_next_parameters();
        model_collection = model_collection.learn_optimized_structures();
        model_collection = model_collection.do_parameter_estimation(...
            params.BCFW_max_iterations, params.BCFW_fval_epsilon,...
            params.compute_true_logZ, params.reweight_denominator,...
            params.printInterval, params.printTest, params.MaxTime)
        end
    end
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%SAVE_RESULTS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      if shuffle_status == 0
       [best_model_index] = get_best_model(model_collection);
       [best_model] = SingleLoopyModel(model_collection, best_model_index);
       best_model=struct(best_model);
       model_collection = struct(model_collection);
       save(strcat(params.exptdir, '/', 'model_collection.mat'), 'model_collection');
       save(strcat(params.exptdir, '/', 'best_model.mat'), 'best_model');
      end
      
end
      