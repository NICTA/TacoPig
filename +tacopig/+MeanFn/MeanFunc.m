% Mean Function Abstract Class

classdef MeanFunc < handle
    
    methods(Abstract)
        
        % Automatically compute hyperparameter dimensionality
        n = npar(D);

        % Mean Function Itself
        mu = eval(X, par);
    end
    
    methods

        % Gradient stub that may be overloaded
        function gradient(this) 
            error([class(this),' does not implement gradients!']);
        end
    end
    
end    