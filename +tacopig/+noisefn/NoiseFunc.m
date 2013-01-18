% Noise Function Abstract Class
% All noise function classes must inherent from this class.

classdef NoiseFunc < tacopig.taco

    methods (Abstract)
        % Returns the number of parameters required by the mean function.
        % This is an abstract method.
        % n_theta = npar(this, D); 
        % Input :   D (dimensionality of the input data)
        % Output :  n (the number of parameters required by the mean function)       
        n_theta = npar(this, D); 
        % Evaluate the noise function at location X given parameters theta
        % This is an abstract method
        %   noise = eval(this, X, theta);
        % Inputs : X (input locations) D x N
        %          theta (parameters) 1 x Number of parameters
        % Outpus : noise (the noise at the input locations) 1 x N
        noise = eval(this, X1, theta);
        
    end
    
    methods (Static, Hidden = true)
        
         function theta = getNoisePar(GP) % returns the number of noise parameters
            if isa(GP, 'tacopig.gp.GpCore')
                theta = GP.noisepar; % Included to help with auto-testing
            elseif isa(GP, 'double')
                theta = GP;
            else
                error('tacopig:badConfiguration', 'Error interpreting noisepar.');
            end
         end
        
    end
    
    methods
        function gradient(this, varargin)
          % Returns the noise gradient with respect to the parameters. Stub that may be overloaded
            error([class(this),' does not implement gradients!']);
        end
        
       
        
    end
end    
