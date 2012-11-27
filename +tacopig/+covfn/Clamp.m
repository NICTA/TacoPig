% Allows a projection matrix like G to be applied around the kernel
classdef GP_ClampCov < GP_CovFunc
    
    properties
       indx
       value
       covfn
    end
    
    methods
        
        function this = GP_ClampCov(covfn, indx, value)
           if ~isa(covfn,'GP_CovFunc')
               error([class(this),': must specify a valid covariance function']); 
           end
           
           this.indx = indx;
           this.value = value;
           this.covfn = covfn;
        end   
            
        function n_theta = npar(this,D)
            n_theta = this.covfn.npar(D) - length(this.indx);
        end
        
        function K = eval(this, X1, X2, parin)
            par = this.getpar(parin, size(X1,1));
            K = this.covfn.eval(X1,X2,par);
        end
        
        function K = Keval(this, X, parin)
            par = this.getpar(parin, size(X,1));
            K = this.covfn.Keval(X,par);
        end
        
        function g = gradient(this,X, parin)
            % Inefficiency at the cost of being clean
            % The alternative is to have a gradient method that takes as an
            % input a list of which gradients we need to know...
            
            [par, nindx] = this.getpar(parin, size(X,1));
            g0 = this.covfn.gradient(X, par); % be careful to concatenate along the right dimension...
            g = g0(nindx);
            
        end
        
        % Overload the point covariance - its trivial to add them
        function v = pointval(this, x_star, parin)
            par = this.getpar(parin, size(x_star,1));
            v = this.covfn.pointval(x_star, par);
        end
        
        function [par,nindx] = getpar(this, parin, D)
            indx = this.indx;
            npar = this.covfn.npar(D);
            par = ones(1,npar);
            par(indx) = 0;
            nindx = find(par>0); % not needed here but used in gradient
            par(nindx) = parin;
            par(indx) = this.value;
        end
        
    end
end
