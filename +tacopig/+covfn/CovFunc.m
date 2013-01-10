% Covariance Function Abstract Class

classdef CovFunc < tacopig.taco
    
    methods (Abstract)
        % Automatically report hyperparameter dimensionality
        % 
        % n_theta = GP.covfn.npar(D); 
        %
        % Gp.CovFn is an instantiation of the covariance function class
        % Inputs: D = dimensionality of the dataset
        % Outputs: n_theta = number of hyperparameters required
        n_theta = npar(this, D); 
        
        % Get covariance matrix between input sets X1,X2
        %
        % K = Gp.CovFn.eval(X1, X2, theta);
        %
        % Gp.CovFn is an instantiation of the covariance function class 
        % X1    - D x N input points 1
        % X2    - D x M input points 2
        % theta - 1 x N vector of hyperparameters
        K = eval(this, X1, X2, theta);
    end

    methods

        % just run with the default constructor
        
        
        function v = pointval(this, x_star, GP)
        % Efficient case for getting diag(k(X,X))
        %
        %         v = GP.covfn.npar.pointval(X, GP)
        %
        % Gp.CovFn is an instantiation of the covariance function class 
        % Inputs : X    - D x N input points            
        %          GP = The GP class instance can be passed to give the covariance function access to its properties
        % Outputs : v = diag(k(X,X)) (The diagnal of the covariance between)
        
            nstar = size(x_star,2);
            v = zeros(1,nstar);
            for i=1:nstar
                xx = x_star(:,i);
                v(i) = this.eval(xx,xx, GP); % using the case where it is a double
            end
        end
   
        
        
        function K = Keval(this, X, GP)
        % Evaluation of k(X,X) (symmetric case)
        %
        % K = GP.covfn.Keval(X, GP)
        % 
        % Inputs :  X = D x N input points
        %           GP = The GP class instance can be passed to give the covariance function access to its properties
        % Outputs : K = covariance matrix between input sets X and itself (N x N)   
            K = this.eval(X,X,GP);
        end
        
        
        function gradient(this)
            % Returns gradient of k(X,X) with respect to each hyperparameter
            error('tacopig:badConfiguration',[class(this),' does not implement gradients!']);
        end
        
        function theta = getCovPar(this, GP)
            % Returns the covariance function's hyperparameters
            if isa(GP, 'tacopig.gp.GpCore')
                theta = GP.covpar;
            elseif isa(GP, 'double')
                theta = GP;
            else
                error('tacopig:badConfiguration', 'Error interpreting covpar.');
            end
        end
        
        
    end
end    
