% Covariance Function Abstract Class

classdef CovFunc < handle
    
    methods (Abstract)
        % Automatically report hyperparameter dimensionality
        n_theta = npar(this, D); 
        
        % Get covariance matrix between input sets X1,X2
        % X1    - D x N input points 1
        % X2    - D x M input points 2
        % theta - 1 x N vector of hyperparameters
        K = eval(this, X1, X2, theta);
    end

    methods
        % Additional Special Cases:        
        
        % Efficient case for getting Kx*x* (just points)
        function v = pointval(this, x_star, theta)
            nstar = size(x_star,2);
            v = zeros(1,nstar);
            for i=1:nstar
                xx = x_star(:,i);
                v(i) = this.eval(xx,xx, theta);
            end
        end
        
        % Evaluation of K(X,X) (symmetric case)
        function K = Keval(this, X, theta)
            K = this.eval(X,X,theta);
        end
        
        % Gradient of K(X,X) (stub)
        function gradient(this)
            error('tacopig:badConfiguration',[class(this),' does not implement gradients!']);
        end
        
    end
end    
