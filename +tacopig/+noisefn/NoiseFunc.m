% Noise Function Abstract Class

classdef NoiseFunc < handle

    methods (Abstract)
        % Automatically report number of parameters:
        n_theta = npar(this); 
        
        % Evaluate
        noise = eval(this, X1, X2, theta);
    end
    
    methods
        
        % Gradient stub that may be overloaded
        function gradient(this)
            error([class(this),' does not implement gradients!']);
        end
    end
end    
