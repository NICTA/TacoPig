% Noise Function Abstract Class

classdef NoiseFunc < tacopig.taco

    methods (Abstract)
        % Automatically report number of parameters:
        n_theta = npar(this, D); 
        % Evaluate
        noise = eval(this, X1, theta);
        
    end
    
    methods (Static)
        
         function theta = getNoisePar(GP)
            if isa(GP, 'tacopig.gp.GpCore')
                theta = GP.noisepar;
            elseif isa(GP, 'double')
                theta = GP;
            else
                error('tacopig:badConfiguration', 'Error interpreting noisepar.');
            end
         end
        
    end
    
    methods
        % Gradient stub that may be overloaded
        function gradient(this, varargin)
            error([class(this),' does not implement gradients!']);
        end
        
       
        
    end
end    
