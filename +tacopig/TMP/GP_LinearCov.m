% Get a covariance function that is well suited to adjacent
% blocks of unknown translation axis-wise

classdef GP_LinearCov < GP_CovFunc
    
    % Most covariance functions will be static
    methods(Static) 
        
        function n_theta = npar(D)
            n_theta = D;       % distance to zero in each dimension
        end
        
        function K = eval(X1, X2, par)
            [D,N1] = size(X1); % number of points in X1
            N2 = size(X2,2);   % number of points in X2
            if D~=size(X2,1)
                error('Dimensionality of X1 and X2 must be the same');
            end
            if (length(par)~=D)
                error('Wrong number of hyperparameters for LinearCov');
            end
            
            % Compute weighted distances:
            w = abs(par);
            K = ones(N1,N2);
            for i=1:D % all input dimensions
                xi = X1(i,:);
                xj = X2(i,:);
                K = K .* max(0,w(i) -abs(xi(ones(1,N2),:)' - xj(ones(1,N1),:)))./w(i);
            end
            
        end
                
        % Also overload the point covariance kx*x* - try this smoother
        % overlap based covariance function
        function v = pointval(x_star, par)
            v = ones(1,size(x_star,2));
        end
        
    end
end
