% Defines a covariance function sum of the children.
% tacopig.covfn.Sum(child1, child2, ...)
% All children must inherit tacopig.covfn.CovFn

classdef Sum < tacopig.covfn.CovFunc
    
    properties(Constant)
        ExampleUsage = 'tacopig.covfn.Sum(tacopig.covfn.Mat3(), tacopig.covfn.SqExp())'; %Instance of class created for testing
    end
    
    properties
       Children % A cell array of the covariance functions that are to be summed
    end
    
    methods
        
        function this = Sum(varargin)
            % Sum covariance function constructor
            % GP.CovFn = tacopig.covfn.Sum(covfunc1, covfunc2, ...)
            % Inputs: covfunc1 (An instance of a covariance function that will be summed to the others

           n_children = length(varargin);
           for i=1:n_children
               if ~isa(varargin{i},'tacopig.covfn.CovFunc')
                  error('Argument ', num2str(i), ' is not a valid covariance function'); 
               end
           end
           this.Children = varargin;
        end   
            
        function n_theta = npar(this,D)
        % Returns the number of hyperparameters required by all of the children covariance functions
        %
        % n_theta = GP.covfunc.npar(D)
        %
        % Gp.CovFn is an instantiation of the Sum covariance function class
        % Inputs: D (Dimensionality of the dataset)
        % Outputs: n_theta (number of hyperparameters required by all of the children covariance function)
        
            n_children = length(this.Children);
            n_theta = 0;
            for i=1:n_children
                n_theta = n_theta + this.Children{i}.npar(D);
            end
        end
        
        function K = eval(this, X1, X2, GP)
        %Get covariance matrix between input sets X1,X2 
        %
        % K = GP.covfunc.eval(X1, X2, GP)
        %
        % Gp.CovFn is an instantiation of the Sum covariance function class
        % Inputs:  X1 = Input locations (D x N matrix)
        %          X2 = Input locations (D x M matrix)
        %          GP = The GP class instance can be passed to give the covariance function access to its properties
        % Outputs: K = covariance matrix between input sets X1,X2 (N x M)
        
            par = this.getCovPar(GP);
            D = size(X1,1); %number of points in X1
            if D~=size(X2,1)
                 error('tacopig:dimMismatch','Dimensionality of X1 and X2 must be the same.');
            end
            
            
            npar = length(par);
            n_children = length(this.Children);
            
            
            if (npar~=this.npar(D))
                error('tacopig:inputInvalidLength','Wrong number of hyperparameters for NegExp');
            end
            
            
            
            K = 0;
            left = 0;
            right = 0;
            for i=1:n_children
                right = right + this.Children{i}.npar(D);
                if (npar<right)
                    error('Need more hyperparameters for CovSum');
                end
                K = K + this.Children{i}.eval(X1, X2, par(left+1:right));
                left = right;
            end
            
        end
        
        
         function K = Keval(this, X, GP)
         % Evaluation of k(X,X) (symmetric case)
         % 
         % K = Keval(X, GP)
         %
         % Gp.CovFn is an instantiation of the Sum covariance function class
         %
         % Inputs:  X = Input locations (D x N matrix)
         %          GP = The GP class instance can be passed to give the covariance function access to its properties
         % Outputs: K = covariance matrix between input sets X and itself (N x N)
            
            par = this.getCovPar(GP);
            D = size(X,1); %number of points in X1
            npar = length(par);
            if (npar~=this.npar(D))
                error('tacopig:inputInvalidLength','Wrong number of hyperparameters for NegExp');
            end
            
            n_children = length(this.Children);
            K = 0;
            left = 0;
            right = 0;
            for i=1:n_children
                right = right + this.Children{i}.npar(D);
                if (npar<right)
                    error('Need more hyperparameters for CovSum');
                end
                K = K + this.Children{i}.Keval(X, par(left+1:right));
                left = right;
            end
            
        end
        
        function g = gradient(this,X, GP)
         % Returns gradient of the covariance matrix k(X,X) with respect to each hyperparameter
        % g = Gp.CovFn.gradient(X, GP)
        %
        % Gp.CovFn is an instantiation of the Sum covariance function class
        %
        % Inputs:   X = Input locations (D x N matrix)
        %           GP = The GP class instance can be passed to give the covariance function access to its properties
        % Outputs:  g = gradient of the covariance matrix k(X,X) with respect to each hyperparameter. Cell Array (1 x Number of Hyperparameters). Each Element is a N x N matrix

            par = this.getCovPar(GP);
            [D,N] = size(X);
            npar = length(par);
            n_children = length(this.Children);
            K = 0;
            left = 0;
            right = 0;
            glist = cell(n_children,1);
            for i=1:n_children
                right = right + this.Children{i}.npar(D);
                if (npar<right)
                    error('Need more hyperparameters for CovSum');
                end
                glist{i} = this.Children{i}.gradient(X, par(left+1:right) );
                left = right;
            end
            g = cat(2,glist{:}); % be careful to concatenate along the right dimension...
            
        end
        
        % Overload the point covariance - its trivial to add them
        function v = pointval(this, x_star, GP)
            % Efficient case for getting diag(k(X,X))
            % v = Gp.CovFn.pointval(X, GP)
            %
            % Gp.CovFn is an instantiation of the Sum covariance function class
            %
            % Inputs : X = Input locations (D x n matrix)
            % Output : v = diag(k(X,X))
            
            par = this.getCovPar(GP);
            [D,N] = size(x_star);
            npar = length(par);
            if (npar~=this.npar(D))
                error('tacopig:inputInvalidLength','Wrong number of hyperparameters for NegExp');
            end
            v = 0;
            n_children = length(this.Children);
            K = 0;
            left = 0;
            right = 0;
            for i=1:n_children
                right = right + this.Children{i}.npar(D);
                if (npar<right)
                    error('Need more hyperparameters for CovSum');
                end
                v = v + this.Children{i}.pointval(x_star, par(left+1:right));
                left = right;
            end
        end
        
    end
end
