% Mean Function Abstract Class

classdef MeanFunc < tacopig.taco
    
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
        
        function theta = getMeanPar(this, GP)
            if isa(GP, 'tacopig.gp.GpCore')
                theta = GP.meanpar;
            elseif isa(GP, 'double')
                theta = GP;
            else
                error('tacopig:badConfiguration', 'Error interpreting meanpar.');
            end
         end
        
    end
    
    
    
    
end    