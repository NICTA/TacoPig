% tacopig.covfn.Product(child1, child2, ...)
% Defines a covariance function product of the children.
% All children must inherit tacopig.covfn.CovFn

classdef Product < tacopig.covfn.CovFunc
    
    properties(Constant)
        ExampleUsage = 'tacopig.covfn.Product(tacopig.covfn.SqExp(), tacopig.covfn.Mat3())'; %Instance of class created for testing
    end
    
    properties
       Children     % A cell array of the covariance functions that are to be multiplied
    end
    
    methods
        function this = Product(varargin) 
        % Product covariance function constructor
           n_children = length(varargin);

           for i=1:n_children
               if ~isa(varargin{i}, 'tacopig.covfn.CovFunc')
                  error('tacopig:badConfiguration',[ 'Argument ', num2str(i), ' is not a valid covariance function']); 
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
            D = size(X1,1); % number of points in X1
            if D~=size(X2,1)
                error('tacopig:dimMismatch', 'Dimensionality of X1 and X2 must be the same');
            end
            
            if (size(par,2)~= this.npar(D))
                error('tacopig:inputInvalidLength','Wrong Hyperparameter Length!');
            end
            
            npar = length(par);
            n_children = length(this.Children);
            K = 1;
            left = 0;
            right = 0;
            for i=1:n_children
                right = right + this.Children{i}.npar(D);
                if (npar<right)
                    error('tacopig:inputInvalidLength', 'Need more hyperparameters for product.');
                end
                K = K .* this.Children{i}.eval(X1, X2, par(left+1:right));
                left = right;
            end
            if (right<npar)
                error('tacopig:inputInvalidLength', 'Extra hyperparameters given to product?');
            end
            
        end
        
        
         function K = Keval(this, X, GP)
            % Evaluation of k(X,X) (symmetric case)
            par = this.getCovPar(GP);
            D = size(X,1); %number of points in X1
            npar = length(par);
            n_children = length(this.Children);
            K = 1;
            left = 0;
            right = 0;
            for i=1:n_children
                right = right + this.Children{i}.npar(D);
                if (npar<right)
                    error('tacopig:inputInvalidLength', 'Need more hyperparameters for product.');
                end
                K = K .* this.Children{i}.Keval(X, par(left+1:right));
                left = right;
            end
            if (right<npar)
                error('tacopig:inputInvalidLength', 'Extra hyperparameters given to product?');
            end
            
        end
        
        function g = gradient(this,X, GP)
            % Returns gradient of k(X,X) with respect to each hyperparameter
            par = this.getCovPar(GP);
            [D,N] = size(X);
            npar = length(par);
            n_children = length(this.Children);
            K = 1;
            left = 0;
            right = 0;
            glist = cell(n_children,1);
            for i=1:n_children
                right = right + this.Children{i}.npar(D);
                if (npar<right)
                    error('tacopig:inputInvalidLength', 'Need more hyperparameters for product.');
                end
                glist{i} = this.Children{i}.gradient(X, par(left+1:right) );
                klist{i} = this.Children{i}.Keval(X, par(left+1:right) );
                left = right;
            end
            
            % Evaluates gradient using:
            % d(ABC)/dtheta = (dA/dtheta * BC) + (dB/dtheta * AC) + (dC/dtheta * AB) 
            counter = 1;
            for i=1:n_children
                for j=1:this.Children{i}.npar(D)
                    ind = 1:numel(klist);
                    others = setxor(ind,i);
                    kothers = ones(N);
                    for t = 1:numel(others)
                        kothers = kothers.*klist{others(t)};
                    end
                     g{counter} = cell2mat(glist{i,1}(j)).*kothers;
                   counter = counter+1;
                end
            end  
            
        end
        
        % Overload the point covariance - its trivial to add them
        function v = pointval(this, x_star, GP)
            par = this.getCovPar(GP);
            [D,N] = size(x_star);
            npar = length(par);
            v = 1;
            n_children = length(this.Children);
            left = 0;
            right = 0;
            for i=1:n_children
                right = right + this.Children{i}.npar(D);
                if (npar<right)
                    error('tacopig:inputInvalidLength', 'Need more hyperparameters for product');
                end
                v = v .* this.Children{i}.pointval(x_star, par(left+1:right));
                left = right;
            end
            if (right<npar)
                error('tacopig:inputInvalidLength', 'Extra hyperparameters given to product?');
            end
            
        end
        
    end
end
