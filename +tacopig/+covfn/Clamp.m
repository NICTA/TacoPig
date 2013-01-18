% Clamps specific parameters of a covariance function
% These parameters are not altered during the learning phase.
% Usage: tacopig.covfn.Clamp(covfn,indx,value);
% 
% Example instantiation:
% GP.CovFn   = tacopig.covfn.Clamp(tacopig.covfn.SqExp(),[1 3],[4 0.6]);
% Instantiates a Squared exponential covariance function (as a property of a Gaussian process instantiation called GP) 
% with its first and third hyperparameters clamped to the value 4 and 0.6, respectively.


classdef Clamp < tacopig.covfn.CovFunc
    
    properties(Constant)
        ExampleUsage = 'tacopig.covfn.Clamp(tacopig.covfn.SqExp(), 1, 2)'; %Instance of class created for testing
    end
    
    properties
       indx     % index of hypermeter(s) to be clamped
       value    % value(s) of clamped hypermeter(s)
       covfn    % the constructor of the covariance function whose hyperparameter(s) are to be clamped
    end
    
    methods
        
        function this = Clamp(covfn, indx, value)
           % Clamp covariance function constructor
           %
           % GP.CovFn = tacopig.covfn.Clamp(covfn, indx, value)
           %
           % Inputs:   covfn     the constructor of the covariance function whose hyperparameter(s) are to be clamped
           %           indx      index of hypermeter(s) to be clamped
           %           value     value(s) of clamped hypermeter(s)
       
           if ~isa(covfn,'tacopig.covfn.CovFunc')
               error('tacopig:badConfiguration', [class(this),': must specify a valid covariance function']); 
           end
           this.indx = indx;
           this.value = value;
           this.covfn = covfn;
        end   
            
        function n_theta = npar(this,D)
            % Returns the number of free hyperparameters required by the covariance function
            %
            % n_theta = GP.CovFn.npar(D)
            %
            % GP.CovFn is an instantiation of the Clamp covariance function class
            % Inputs: D = Dimensionality of the input data
            % Outputs: n_theta = the number of free hyperparameters required by the covariance function
            
            n_theta = this.covfn.npar(D) - length(this.indx);
        end
        
        function K = eval(this, X1, X2, GP)
            %Get covariance matrix between input sets X1,X2
            %
            % K = GP.CovFn.eval(X1, X2, GP)
            %
            % GP.CovFn is an instantiation of the Clamp covariance function class
            % Inputs:  X1 = Input locations (D x N matrix)
            %          X2 = Input locations (D x M matrix)
            %          GP = The GP class instance can be passed to give the covariance function access to its properties
            % Outputs: K = covariance matrix between input sets X1,X2 (N x M)
            
            parin = this.getCovPar(GP);
            if (length(parin)~=this.npar(size(X1,1)))
                error('tacopig:inputInvalidLength','Wrong number of hyperparameters!');
            end
            par = this.getpar(parin, size(X1,1));
            K = this.covfn.eval(X1,X2,par);
        end
        
        function K = Keval(this, X, GP)
            % Evaluation of k(X,X) (symmetric case)
            % Get covariance matrix between input sets X1,X2
            %
            % K = GP.CovFn.Keval(X, GP)
            %
            % GP.CovFn is an instantiation of the Clamp covariance function class
            %
            % Inputs:  X = Input locations (D x N matrix)
            %          GP = The GP class instance can be passed to give the covariance function access to its properties
            % Outputs: K = covariance matrix between input sets X and itself (N x N)            
            parin = this.getCovPar(GP);
            
            if (length(parin)~=this.npar(size(X,1)))
                error('tacopig:inputInvalidLength','Wrong number of hyperparameters!');
            end
            par = this.getpar(parin, size(X,1));
            K = this.covfn.Keval(X,par);
        end
        
        function g = gradient(this,X, GP)
            % Returns gradient of the covariance matrix k(X,X) with respect to each hyperparameter
            % g = GP.CovFn.gradient(X, GP)
            %
            % GP.CovFn is an instantiation of the Clamp covariance function class
            %
            % Inputs:   X = Input locations (D x N matrix)
            %           GP = The GP class instance can be passed to give the covariance function access to its properties
            % Outputs:  g = gradient of the covariance matrix k(X,X) with respect to each hyperparameter. Cell Array (1 x Number of Hyperparameters). Each Element is a N x N matrix


            
            parin = this.getCovPar(GP);
            % Inefficiency at the cost of being clean
            % The alternative is to have a gradient method that takes as an
            % input a list of which gradients we need to know...
            
            [par, nindx] = this.getpar(parin, size(X,1));
            g0 = this.covfn.gradient(X, par); % be careful to concatenate along the right dimension...
            g = g0(nindx);
            
        end
        
        % Overload the point covariance - its trivial to add them
        function v = pointval(this, x_star, GP)
            % Efficient case for getting diag(k(X,X))
            % v = GP.CovFn.pointval(X, GP)
            %
            % GP.CovFn is an instantiation of the Clamp covariance function class
            %
            % Inputs : X = Input locations (D x n matrix)
            %          GP = The GP class instance can be passed to give the covariance function access to its properties
            % Output : v = diag(k(X,X))
            parin = this.getCovPar(GP);
            if (length(parin)~=this.npar(size(x_star,1)))
                error('tacopig:inputInvalidLength','Wrong number of hyperparameters!');
            end
            par = this.getpar(parin, size(x_star,1));
            v = this.covfn.pointval(x_star, par);
        end
    end
    
    methods(Access = private)
        
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
