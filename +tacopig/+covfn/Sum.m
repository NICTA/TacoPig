% Defines a covariance function sum of the children.
% tacopig.covfn.Sum(child1, child2, ...)
% All children must inherit tacopig.covfn.CovFn

classdef Sum < tacopig.covfn.CovFunc
    
    properties(Constant)
        ExampleUsage = 'Sum(tacopig.covfn.Mat3(), tacopig.covfn.SqExp())'; %Instance of class created for testing
    end
    
    properties
       Children % A cell array of the covariance functions that are to be summed
    end
    
    methods
        
        function this = Sum(varargin)
        % Sum covariance function constructor
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
            n_children = length(this.Children);
            n_theta = 0;
            for i=1:n_children
                n_theta = n_theta + this.Children{i}.npar(D);
            end
        end
        
        function K = eval(this, X1, X2, GP)
        %Get covariance matrix between input sets X1,X2 
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
        % Returns gradient of k(X,X) with respect to each hyperparameter
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
