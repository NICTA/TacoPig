% tacopig.covfn.Product(child1, child2, ...)
% Defines a covariance function product of the children.
% All children must inherit tacopig.covfn.CovFn

classdef Product < tacopig.covfn.CovFunc
    
    properties(Constant)
        teststring = 'Product(tacopig.covfn.SqExp(), tacopig.covfn.Mat3())'
    end
    
    properties
       children
    end
    
    methods
        function this = Product(varargin)
           n_children = length(varargin);

           for i=1:n_children
               if ~isa(varargin{i}, 'tacopig.covfn.CovFunc')
                  error('tacopig:badConfiguration',[ 'Argument ', num2str(i), ' is not a valid covariance function']); 
               end
           end
           this.children = varargin;
        end   
            
        function n_theta = npar(this,D)
            n_children = length(this.children);
            n_theta = 0;
            for i=1:n_children
                n_theta = n_theta + this.children{i}.npar(D);
            end
        end
        
        function K = eval(this, X1, X2, GP)
            par = this.getCovPar(GP);
            D = size(X1,1); % number of points in X1
            if D~=size(X2,1)
                error('tacopig:dimMismatch', 'Dimensionality of X1 and X2 must be the same');
            end
            
            if (size(par,2)~= this.npar(D))
                error('tacopig:inputInvalidLength','Wrong Hyperparameter Length!');
            end
            
            npar = length(par);
            n_children = length(this.children);
            K = 1;
            left = 0;
            right = 0;
            for i=1:n_children
                right = right + this.children{i}.npar(D);
                if (npar<right)
                    error('tacopig:inputInvalidLength', 'Need more hyperparameters for product.');
                end
                K = K .* this.children{i}.eval(X1, X2, par(left+1:right));
                left = right;
            end
            if (right<npar)
                error('tacopig:inputInvalidLength', 'Extra hyperparameters given to product?');
            end
            
        end
        
        
         function K = Keval(this, X, GP)
            par = this.getCovPar(GP);
            D = size(X,1); %number of points in X1
            npar = length(par);
            n_children = length(this.children);
            K = 1;
            left = 0;
            right = 0;
            for i=1:n_children
                right = right + this.children{i}.npar(D);
                if (npar<right)
                    error('tacopig:inputInvalidLength', 'Need more hyperparameters for product.');
                end
                K = K .* this.children{i}.Keval(X, par(left+1:right));
                left = right;
            end
            if (right<npar)
                error('tacopig:inputInvalidLength', 'Extra hyperparameters given to product?');
            end
            
        end
        
        function g = gradient(this,X, GP)
            par = this.getCovPar(GP);
            [D,N] = size(X);
            npar = length(par);
            n_children = length(this.children);
            K = 1;
            left = 0;
            right = 0;
            glist = cell(n_children,1);
            for i=1:n_children
                right = right + this.children{i}.npar(D);
                if (npar<right)
                    error('tacopig:inputInvalidLength', 'Need more hyperparameters for product.');
                end
                glist{i} = this.children{i}.gradient(X, par(left+1:right) );
                klist{i} = this.children{i}.Keval(X, par(left+1:right) );
                left = right;
            end
            
            % Evaluates gradient using:
            % d(ABC)/dtheta = (dA/dtheta * BC) + (dB/dtheta * AC) + (dC/dtheta * AB) 
            counter = 1;
            for i=1:n_children
                for j=1:this.children{i}.npar(D)
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
            n_children = length(this.children);
            left = 0;
            right = 0;
            for i=1:n_children
                right = right + this.children{i}.npar(D);
                if (npar<right)
                    error('tacopig:inputInvalidLength', 'Need more hyperparameters for product');
                end
                v = v .* this.children{i}.pointval(x_star, par(left+1:right));
                left = right;
            end
            if (right<npar)
                error('tacopig:inputInvalidLength', 'Extra hyperparameters given to product?');
            end
            
        end
        
    end
end
