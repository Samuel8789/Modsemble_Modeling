classdef SingleLoopyModel
    %SINGLELOOPYMODEL On single model (generated by one combination of
    %regularization parameters and density).
    
    properties
        % training and test sets
        x_train, x_test;
        
        % struct with 
        %       F (2xnode_count matrix with node weights) 
        %       G (4xedge_count matrix with edge weights), 
        %       node_potentials (column vector),
        %       edge_potentials (symmetric matrix),
        %       true_logZ (computed with JTA)
        %       logZ (computed with Bethe approx.)
       theta;
        
        %   s_lambda: regularization of structure learning
        %   p_lambda: regularization of parameter learning
        %   density: target density of the structure
        s_lambda, p_lambda, density;
        
        %   structure: NxN binary adjacency matrix
        structure,
        
        %   max_degree: maximum number of connection to a node
        max_degree, median_degree, mean_degree, rms_degree;
        reweight;
        
        %   train_likelihood: avg. (per sample) likelihood of training set
        %   test_likelihood: avg. (per sample) likelihood of test set
        train_likelihood,
        test_likelihood,
        
        ep_on = [];

    end
    
    methods
        
        % Constructor: used by LoopyModelCollection objects
        function self = SingleLoopyModel(model_collection, best_model_index)
            self.x_train = model_collection.x_train;
            self.x_test = model_collection.x_test;
            model_struct = model_collection.models{best_model_index};
            self.theta = model_struct.theta;
            self.s_lambda = model_struct.s_lambda;
            self.p_lambda = model_struct.p_lambda;
            self.structure = model_struct.structure;
            self.train_likelihood = model_struct.train_likelihood;
            self.test_likelihood = model_struct.test_likelihood;
            self.max_degree = model_struct.max_degree;
            
            if isfield(model_struct, 'median_degree')
                self.median_degree = model_struct.median_degree;
            end
            if isfield(model_struct, 'mean_degree')
                self.mean_degree = model_struct.mean_degree;
            end
            if isfield(model_struct, 'rms_degree')
                self.rms_degree = model_struct.rms_degree;
            end
            if isfield(model_struct, 'reweight')
                self.reweight = model_struct.reweight;
            end
            
        end
    end
end    
  

