% Gaussian Process Core Model
% Abstract class to define basic GP interface.

classdef GpCore < tacopig.taco
    
    % Members common to all GP models
    properties
        MeanFn          % An instance of a mean function class
        CovFn           % An instance of a covariance function class
        NoiseFn         % An instance of a noise function class
        X               % Training inputs
        y               % Target vector
    end
    
    % Methods that all GP models are required to implement
    methods (Abstract)
        
        % Automatic model selection. (Required method for all GP classes)
        learn(this);                        
        
        % Probabilistic inference. (Required method for all GP classes)
        solve(this);
        
        % Batch querying. (Required method for all GP classes)
        [mu, var] = query(this, xstar);
        
        % Objective function for learning. (Required method for all GP classes)
        [o, ograd] = objectfun(theta);
        
    end
end    