% Constant Mean GP Function
% This is an example implementation of a polynomial mean function...

classdef StationaryMean < tacopig.meanfn.MeanFunc
    
    
    methods(Static)
        function n = npar(~) 
            n = 1; % intercept and gradients in each dimension
        end
        
        function mu = eval(X, par) 
            if (numel(par)~=1)
                error('tacopig:inputInvalidLength','wrong number of hyperparameters!')
            end
            
            N = size(X,2);
            mu = par(1)*ones(1,N);
        end
        
        function g = gradient(X, ~) 
            N = size(X,2);
            g = cell(1,1);
            g{1} = ones(1,N);
        end
    end
end    


