% Mean Function Abstract Class
% All mean function classes must inherent from this class.

classdef MeanFunc < tacopig.taco
    
    methods(Abstract)
        
        % Returns the number of parameters required by the mean function.
        % n = npar(D); 
        % Input :   D (dimensionality of the input data)
        % Output :  n (the number of parameters required by the mean function)
        n = npar(D);         


        % Returns the value of the mean function.
        % mu = eval(X, par); 
        % Input : X (input points D x N)
        %         par (mean function parameters).
        % Output : mu (the mean function at the location of the input points)
        mu = eval(X, par);   


    end
    
    methods

        function gradient(this) 
        % Returns the mean gradient with respect to the parameters. Stub that may be overloaded
            error([class(this),' does not implement gradients!']);
        end
        
        function theta = getMeanPar(this, GP)
        % Returns the parameters of the mean function    
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