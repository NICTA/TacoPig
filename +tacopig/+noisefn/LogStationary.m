classdef LogStationary < tacopig.noisefn.NoiseFunc
   
    methods(Static) 
        
        function n_theta = npar(~)
            n_theta = 1; % normal iid noise
        end
        
        function noise = eval(X, GP)
            par = tacopig.noisefn.NoiseFunc.getNoisePar(GP);
            [D,N] = size(X); %N number of points in X
            noise = exp(par)*eye(N);
        end
        
        function [g] = gradient(X, GP)
            par = tacopig.noisefn.NoiseFunc.getNoisePar(GP);
            [D,N] = size(X); %N number of points in X
            g = {exp(par)*eye(N)};
        end

    end
end
