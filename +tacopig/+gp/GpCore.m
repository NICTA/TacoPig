% Gaussian Process Core Model
% Abstract class to define basic GP interface.

classdef GpCore < tacopig.taco
    
    % Members common to all GP models
    properties
        MeanFn
        CovFn
        NoiseFn
        X
        y
    end
    
    % Methods that all GP models are required to implement
    methods (Abstract)
        learn(this);                        % Automatic model selection
        solve(this);                        % Probabilistic inference
        [mu, var] = query(this, xstar);     % Batch querying
        [o, ograd] = objectfun(theta);      % Objective function for learning
    end
end    