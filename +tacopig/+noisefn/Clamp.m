% Allows a clamp hyperparameters in the noise func
classdef Clamp < tacopig.noisefn.NoiseFunc
    
    properties
       indx
       value
       noisefn
    end
    
    properties(Constant)
        ExampleUsage = 'tacopig.noisefn.Clamp(tacopig.noisefn.Stationary(), 1, 0)';
    end
    
    methods
        
        function this = Clamp(noisefn, indx, value)
           if ~isa(noisefn,'tacopig.noisefn.NoiseFunc')
               error([class(this),': must specify a valid noise function']); 
           end
           
           this.indx = indx;
           this.value = value;
           this.noisefn = noisefn;
        end   
            
        function n_theta = npar(this,D)
            n_theta = this.noisefn.npar(D) - length(this.indx);
        end
        
        function noise = eval(this, X, GP)
            parin = this.getNoisePar(GP);
            par = this.getpar(parin);
            noise = this.noisefn.eval(X,par);
        end
        
        
        function g = gradient(this,X, GP)
            parin = this.getNoisePar(GP);
            [par, nindx] = this.getpar(parin);
            g0 = this.noisefn.gradient(X, par); % be careful to concatenate along the right dimension...
            g = g0(nindx);
        end
        
 
        function [par,nindx] = getpar(this, parin)
            indx = this.indx;
            npar = this.noisefn.npar;
            par = ones(1,npar);
            par(indx) = 0;
            nindx = find(par>0); % not needed here but used in gradient
            par(nindx) = parin;
            par(indx) = this.value;
        end
        
    end
end
