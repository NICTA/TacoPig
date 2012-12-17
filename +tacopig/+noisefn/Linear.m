classdef Linear < tacopig.noisefn.NoiseFunc
   
    % Most noise functions will be static
    methods(Static) 
        
        function n_theta = npar(~)
            n_theta = 1;
        end
        
        function noise = eval(X, par)
             noise = diag(abs(X(1,:)*(par)^2));
        end
        
        function [g] = gradient(X, par)
             g = {2*diag(abs(X(1,:)))*par};
        end
    end
    
end
