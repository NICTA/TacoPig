% Clamps specific parameters of a covariance function so that they are not
% altered during the learning phase.

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
           if ~isa(covfn,'tacopig.covfn.CovFunc')
               error('tacopig:badConfiguration', [class(this),': must specify a valid covariance function']); 
           end
           this.indx = indx;
           this.value = value;
           this.covfn = covfn;
        end   
            
        function n_theta = npar(this,D)
            % Returns the number of free hyperparameters required by the covariance function
            n_theta = this.covfn.npar(D) - length(this.indx);
        end
        
        function K = eval(this, X1, X2, GP)
            %Get covariance matrix between input sets X1,X2
            parin = this.getCovPar(GP);
            if (length(parin)~=this.npar(size(X1,1)))
                error('tacopig:inputInvalidLength','Wrong number of hyperparameters!');
            end
            par = this.getpar(parin, size(X1,1));
            K = this.covfn.eval(X1,X2,par);
        end
        
        function K = Keval(this, X, GP)
            % Evaluation of k(X,X) (symmetric case)
            parin = this.getCovPar(GP);
            if (length(parin)~=this.npar(size(X,1)))
                error('tacopig:inputInvalidLength','Wrong number of hyperparameters!');
            end
            par = this.getpar(parin, size(X,1));
            K = this.covfn.Keval(X,par);
        end
        
        function g = gradient(this,X, GP)
            % Returns gradient of k(X,X) with respect to each hyperparameter
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
