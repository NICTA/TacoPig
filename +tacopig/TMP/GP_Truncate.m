% Warning: This covariance function is not positive definite - do not use
% for a normal GP implementation - its purpose is to test the behaviour of 
% SVD based GP learning when we approximate the full covariance function
% with one that is sparse non-psd.

classdef GP_Truncate < GP_CovFunc

    
    properties
        covfn
        threshold
        indx
    end
        % Most covariance functions will be static
    methods
        
        function this = GP_Truncate(covfn, thresh,indx)
            this.covfn = covfn;
            this.threshold = thresh; 
            this.indx = indx;
        end
        
        function n_theta = npar(this,D)
            n_theta = this.covfn.npar(D);
        end
        
        function K = eval(this, X1, X2, par)
            thr = this.threshold*(par(this.indx)^2);
            K = this.covfn.eval(X1, X2, par);
            K(abs(K)<thr) = 0;
        end
        
        function [g] = gradient(this,X, par)
            thr = this.threshold*(par(this.indx)^2);
            K = this.covfn.Keval(X, par);
            g = this.covfn.gradient(X, par);
            setzro = abs(K)<thr;
            for i=1:length(g)
                g{i}(setzro) = 0;
            end
        end
        
        % Also overload the point covariance kx*x* - its trivial
        function v = pointval(this, x_star, par)
            thr = this.threshold*(par(this.indx)^2);
            v = this.covfn.pointval(x_star, par);
            v(abs(v)<thr) = 0;
        end
        
    end
end
