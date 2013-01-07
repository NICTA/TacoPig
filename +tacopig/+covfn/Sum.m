% Covariance function to add two or more underlying covariance functions.

classdef Sum < tacopig.covfn.CovFunc
    
    properties(Constant)
        teststring = 'Sum(tacopig.covfn.Mat3(), tacopig.covfn.SqExp())';
    end
    
    properties
       children
    end
    
    methods
        
        function this = Sum(varargin)
           n_children = length(varargin);
           for i=1:n_children
               if ~isa(varargin{i},'tacopig.covfn.CovFunc')
                  error('Argument ', num2str(i), ' is not a valid covariance function'); 
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
            D = size(X1,1); %number of points in X1
            if D~=size(X2,1)
                 error('tacopig:dimMismatch','Dimensionality of X1 and X2 must be the same.');
            end
            
            
            npar = length(par);
            n_children = length(this.children);
            
            
            if (npar~=this.npar(D))
                error('tacopig:inputInvalidLength','Wrong number of hyperparameters for NegExp');
            end
            
            
            
            K = 0;
            left = 0;
            right = 0;
            for i=1:n_children
                right = right + this.children{i}.npar(D);
                if (npar<right)
                    error('Need more hyperparameters for CovSum');
                end
                K = K + this.children{i}.eval(X1, X2, par(left+1:right));
                left = right;
            end
            
        end
        
        
         function K = Keval(this, X, GP)
            par = this.getCovPar(GP);
            D = size(X,1); %number of points in X1
            npar = length(par);
            if (npar~=this.npar(D))
                error('tacopig:inputInvalidLength','Wrong number of hyperparameters for NegExp');
            end
            
            n_children = length(this.children);
            K = 0;
            left = 0;
            right = 0;
            for i=1:n_children
                right = right + this.children{i}.npar(D);
                if (npar<right)
                    error('Need more hyperparameters for CovSum');
                end
                K = K + this.children{i}.Keval(X, par(left+1:right));
                left = right;
            end
            
        end
        
        function g = gradient(this,X, GP)
            par = this.getCovPar(GP);
            [D,N] = size(X);
            npar = length(par);
            n_children = length(this.children);
            K = 0;
            left = 0;
            right = 0;
            glist = cell(n_children,1);
            for i=1:n_children
                right = right + this.children{i}.npar(D);
                if (npar<right)
                    error('Need more hyperparameters for CovSum');
                end
                glist{i} = this.children{i}.gradient(X, par(left+1:right) );
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
            n_children = length(this.children);
            K = 0;
            left = 0;
            right = 0;
            for i=1:n_children
                right = right + this.children{i}.npar(D);
                if (npar<right)
                    error('Need more hyperparameters for CovSum');
                end
                v = v + this.children{i}.pointval(x_star, par(left+1:right));
                left = right;
            end
        end
        
    end
end
