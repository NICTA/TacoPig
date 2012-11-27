% Allows a re-ordering or repetition of hyperparameters
% usage: Remap(covfn, index)
% Eg. [1 1 2] would map optimisation vector [a b] into hyper vector [a a b]

classdef Remap < tacopig.covfn.CovFunc
    
    properties
       indx
       covfn
    end
    
    methods
        
        function this = Remap(covfn, indx)
            
           if ~isa(covfn,'tacopig.covfn.CovFunc')
               error([class(this),': must specify a valid covariance function']); 
           end
           this.indx = indx;
           this.covfn = covfn;
        end   
            
        function n_theta = npar(this,D)
            n_dep = this.covfn.npar(D);
            indx = this.indx;
            if (n_dep~=length(indx))
                error('GP_Remap: index maps to an incorrect number of hypers.');
            end
            if (length(unique(indx)) ~= max(indx))
               error('GP_Remap: must specify an index using consecutive values.')
            end
            n_theta = length(unique(indx));
        end
        
        function K = eval(this, X1, X2, parin)
            par = parin(this.indx);
            K = this.covfn.eval(X1,X2,par);
        end
        
        function K = Keval(this, X, parin)
            par = parin(this.indx);
            K = this.covfn.Keval(X,par);
        end
        
        function g = gradient(this,X, parin)
            
            % Have to be careful and add the gradients here            
            indx = this.indx;
            nindx = length(indx);
            par = parin(indx);
            g0 = this.covfn.gradient(X, par); % be careful to concatenate along the right dimension...
            npar = length(unique(indx));
            g = cell(1, npar);  
            zro = zeros(size(g0{1}));
            for i=1:npar
                g{i} = zro;
            end
            
            for i=1:length(indx) % number of hypers of underlying covariance function.
                ind = indx(i);
                g{ind} = g{ind} + g0{i}; 
            end
            
            
        end
        
        % Overload the point covariance - its trivial to add them
        function v = pointval(this, x_star, parin)
            par = parin(this.indx);
            v = this.covfn.pointval(x_star, par);
        end
        
    end
end
