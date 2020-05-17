
% lasso_node_by_node_group

% ADOPTED AND EDITED BY DARIK ONEIL YUSTE LABORATORY 11/19

% Utilizes glmnet from Stanford Statistics. 

%% Documentation

%function input is:

%x_train: logical, binary matrix i = time and j = node

%variable groups: optional cell array that marks which variable belongs to
%each group to avoid prediction by same group.

%function ouput

%coefficients: matrix where i = lasso and j = coefficient

%% FUNCTION

function [GLM_array] = lasso_node_by_node_group(x_train, variable_groups,LASSO_options)


    options = LASSO_options;
    node_count = size(x_train,2);
   
    disp('Learning Structure')
    GLM_array = cell(node_count,1);
    
   for label_node = 1:node_count
        feature_nodes = variable_groups{label_node};
        X = x_train(:,feature_nodes);
        Y = x_train(:,label_node);
        CVerr = cvglmnet(X,Y,'binomial',options);
        GLM_array{label_node}= CVerr;
        fprintf('.')
   end
  
end

   
