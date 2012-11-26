% Gaussian Process Model Abstract Class

classdef GP_Model < handle
    
    % Member variables that all GP models will have:
    properties
        MeanFn
        CovFn
        NoiseFn
        X
        y
    end
    
    % Methods that all GP models are required to implement
    methods (Abstract)
        theta = learn(this);     % Automatic model selection
        solve(this);             % Probabilistic inference
        [mu, var] = query(this, xstar);  % Batch querying
        
        % Objective for model selection
        % Interestingly, once we change this function pointer it behaves as
        % a member function.
        [objective, objective_grad] = objectfun(theta); 
    end
    
end    