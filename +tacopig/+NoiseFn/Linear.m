classdef GP_LinIncreaseNoise < GP_NoiseFunc
   
    % Most noise functions will be static
    methods(Static) 
        
        function n_theta = npar(~)
            n_theta = 1; % 
        end
        
        function noise = eval(GP, par)
             noise = diag(abs(GP.X(1,:)*(par)^2));
        end
        
        function [g] = gradient(GP, par)
             g = {2*diag(abs(GP.X(1,:)))*par};
        end

    end
end
