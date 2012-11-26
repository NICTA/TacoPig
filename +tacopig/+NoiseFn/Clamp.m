% Allows a clamp hyperparameters in the noise func
classdef GP_ClampNoise < GP_NoiseFunc
    
    properties
       indx
       value
       noisefn
    end
    
    methods
        
        function this = GP_ClampNoise(noisefn, indx, value)
           if ~isa(noisefn,'GP_NoiseFunc')
               error([class(this),': must specify a valid noise function']); 
           end
           
           this.indx = indx;
           this.value = value;
           this.noisefn = noisefn;
        end   
            
        function n_theta = npar(this)
            n_theta = this.noisefn.npar - length(this.indx);
        end
        
        function noise = eval(this, GP, parin)
            par = this.getpar(parin);
            noise = this.noisefn.eval(GP,par);
        end
        
        
        function g = gradient(this,GP, parin)
            % Inefficiency at the cost of being clean
            % The alternative is to have a gradient method that takes as an
            % input a list of which gradients we need to know...
            
            [par, nindx] = this.getpar(parin);
            g0 = this.noisefn.gradient(GP.X, par); % be careful to concatenate along the right dimension...
            
            
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
