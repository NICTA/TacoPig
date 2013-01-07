% Wraps a covariance function as a noise function
% A covariance function isn't accepted by default - you have to explicitly
% request it. We can't wrap all covariance functions
classdef CovNoise < tacopig.noisefn.NoiseFunc
    
    properties(GetAccess = 'public', SetAccess = 'private')
       noisefn
    end
    
    properties(Constant)
        teststring = 'CovNoise(tacopig.covfn.SqExp())';
    end
    
    methods
        
        function this = CovNoise(noisefn)
           if ~isa(noisefn,'tacopig.covfn.CovFunc')
               error([class(this),': wraps a valid covariance function as a noise function!']); 
           end
           this.noisefn = noisefn;
        end   
            
        function n_theta = npar(this,D)
            n_theta = this.noisefn.npar(D);
        end
        
        function noise = eval(this, X, GP)
            parin = this.getNoisePar(GP);
            noise = this.noisefn.Keval(X,parin);
        end
                
        function g = gradient(this,X, GP)
            parin = this.getNoisePar(GP);
            g = this.noisefn.gradient(X, parin);
        end
        
    end
end
