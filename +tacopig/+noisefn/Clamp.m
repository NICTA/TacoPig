% Allows a clamp hyperparameters in the noise func
classdef Clamp < tacopig.noisefn.NoiseFunc
    
    properties
       indx         % index of hypermeter(s) to be clamped
       value        % value(s) of clamped hypermeter(s)
       noisefn      % the constructor of the noisefn function whose hyperparameter(s) are to be clamped
    end
    
    properties(Constant)
        ExampleUsage = 'tacopig.noisefn.Clamp(tacopig.noisefn.Stationary(), 1, 0)'; %Instance of class created for testing
    end
    
    methods
        
        function this = Clamp(noisefn, indx, value)
           % Clamp noise function constructor
           %
           % GP.CovFn = tacopig.noisefn.Clamp(noisefn, indx, value)
           %
           % Inputs:   noisefn   the constructor of the noisefn function whose hyperparameter(s) are to be clamped
           %           indx      index of hypermeter(s) to be clamped
           %           value     value(s) of clamped hypermeter(s)
           
           if ~isa(noisefn,'tacopig.noisefn.NoiseFunc')
               error([class(this),': must specify a valid noise function']); 
           end
           
           this.indx = indx;
           this.value = value;
           this.noisefn = noisefn;
        end   
            
        function n_theta = npar(this,D)
            % Returns the number of free parameters required by the noise function
            %
            % n_theta = GP.NoiseFn.npar(D)
            %
            % GP.NoiseFn is an instantiation of the Clamp noise function class
            % Inputs: D = Dimensionality of the input data
            % Outputs: n_theta = the number of free parameters required by the noise function
           
            n_theta = this.noisefn.npar(D) - length(this.indx);
        end
        
        function noise = eval(this, X, GP)
            % Get noise matrix between the points in input set X
            %
            % noise = GP.NoiseFn.eval(X, GP)
            %
            % GP.CovFn is an instantiation of the Clamp covariance function class
            % Inputs:  X = Input locations (D x N matrix)
            %          GP = The GP class instance can be passed to give the noise function access to its properties
            % Outputs: noise = noise matrix between points in the input set X (N x N)
            
            parin = this.getNoisePar(GP);
            par = this.getpar(parin);
            noise = this.noisefn.eval(X,par);
        end
        
        
        function g = gradient(this,X, GP)
            % Get gradient of the noise function with respect to the parameters
            %
            % g = GP.NoiseFn.gradient(X, GP)
            %
            % GP.NoiseFn is an instantiation of the Clamp covariance function class
            % Inputs:  X = Input locations (D x N matrix)
            %          GP = The GP class instance can be passed to give the noise function access to its properties
            % Outputs: g = gradient of the noise function with respect to the parameters. Cell array, one N x N matrix for each parameter
            
            parin = this.getNoisePar(GP);
            [par, nindx] = this.getpar(parin);
            g0 = this.noisefn.gradient(X, par); % be careful to concatenate along the right dimension...
            g = g0(nindx);
        end
        
    end
    
    methods(Hidden = true)
        function [par,nindx] = getpar(this, parin)
            % Returns the free parameters of the noise function
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
