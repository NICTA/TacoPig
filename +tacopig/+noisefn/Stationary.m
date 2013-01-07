classdef Stationary < tacopig.noisefn.NoiseFunc
   
    % Most noise functions will be static
    methods(Static) 
        
        function n_theta = npar(~)
            n_theta = 1; % normal iid noise
        end
        
        function noise = eval(X, GP)
            par = tacopig.noisefn.NoiseFunc.getNoisePar(GP);
            [D,N] = size(X); %N number of points in X
            noise = par^2*eye(N);
        end
        
        function [g] = gradient(X, GP)
            par = tacopig.noisefn.NoiseFunc.getNoisePar(GP);
            [D,N] = size(X); %N number of points in X
            g = {2*par*eye(N)};
        end

    end
end
