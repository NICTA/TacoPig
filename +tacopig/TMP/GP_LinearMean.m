% Constant Mean GP Function
% This is an example implementation of a polynomial mean function...

classdef GP_LinearMean < GP_MeanFunc

    % Protected member variables:
    properties
    end

    methods(Static)
        function n = npar(D) 
            n = 1+D; % intercept and gradients in each dimension
        end
        
        function mu = eval(X, par) 
            mu = par(1) + (par(2:end)*X);
        end
        
        function g = gradient(X, ~) 
            % neatly our gradient doesnt depend on the parameter values at all
            D = size(X,1);
            N = size(X,2);
            % couldnt g be a matrix?
            g = cell(1,D+1);
            g{1} = ones(1, N);
            for d = 1:D
                g{d+1} = X(d,:);
            end
        end
        
        
    end
end    


