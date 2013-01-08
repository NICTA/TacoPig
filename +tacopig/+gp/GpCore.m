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
        
        % Automatic model selection
        learn(this);                        
        
        % Probabilistic inference
        solve(this);
        
        % Batch querying
        [mu, var] = query(this, xstar);
        
        % Objective function for learning
        [o, ograd] = objectfun(theta);
        
    end
end    