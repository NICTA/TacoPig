% Allows a re-ordering or repetition of hyperparameters
% usage: Remap(covfn, index)
% Eg. [1 1 2] would map optimisation vector [a b] into hyper vector [a a b]

classdef Remap < tacopig.covfn.CovFunc

    properties(Constant)
        ExampleUsage = 'tacopig.covfn.Remap(tacopig.covfn.SqExp(), [1 1 2 3])'; %Instance of class created for testing
    end
    properties
       covfn   % The covariance function to which the remapping of hyperparameters is applied
       indx    % A vector indicating the mapping of the hyperparameter vector to the shorter vector used during learning
    end
    
    methods
        
        function this = Remap(covfn, indx)
        % Constructor  
        %
        % GP.CovFn = tacopig.covfn.Remap(covfn, indx)
        %
        % Inputs:   covfn  = The covariance function to which the remapping of hyperparameters is applied
        %           indx = A vector indicating the mapping of the hyperparameter vector to the shorter vector
        
           if ~isa(covfn,'tacopig.covfn.CovFunc')
               error([class(this),': must specify a valid covariance function']); 
           end
           this.indx = indx;
           this.covfn = covfn;
        end   
            
        function n_theta = npar(this,D)
            % Returns the number of hyperparameters required by the Remap instance
            %
            % n_theta = Gp.CovFn.npar(D)
            %
            % Gp.CovFn is an instantiation of the Remap covariance function class
            %
            % Inputs : D = Dimensionality
            % Outputs: n_theta = number of hyperparameters required by the Remap instance
            
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
        
        function K = eval(this, X1, X2, GP)
            % Get covariance matrix between input sets X1,X2
            %
            % K = Gp.CovFn.eval(X1, X2, GP)
            %
            % Gp.CovFn is an instantiation of the Remap covariance function class
            %
            % Inputs:  X1 = Input locations (D x N matrix)
            %          X2 = Input locations (D x M matrix)
            %          GP = The GP class instance can be passed to give the covariance function access to its properties
            % Outputs: K = covariance matrix between input sets X1,X2 (N x M)
            
            parin = this.getCovPar(GP);
            [D,N1] = size(X1); %number of points in X1
             if (length(parin)~=this.npar(D))
                error('tacopig:inputInvalidLength','Wrong number of hyperparameters!');
            end
            par = parin(this.indx);
            K = this.covfn.eval(X1,X2,par);
        end
        
        function K = Keval(this, X, GP)
            % Evaluation of k(X,X) (symmetric case)
            % Get covariance matrix between input sets X1,X2
            %
            % K = Gp.CovFn.Keval(X, GP)
            %
            % Gp.CovFn is an instantiation of the Remap covariance function class
            %
            % Inputs:  X = Input locations (D x N matrix)
            %          GP = The GP class instance can be passed to give the covariance function access to its properties
            % Outputs: K = covariance matrix between input sets X and itself (N x N)
            parin = this.getCovPar(GP);
            [D,N1] = size(X); %number of points in X1
             if (length(parin)~=this.npar(D))
                error('tacopig:inputInvalidLength','Wrong number of hyperparameters!');
            end
            par = parin(this.indx);
            K = this.covfn.Keval(X,par);
        end
        
        function g = gradient(this,X, GP)
            % Returns gradient of the covariance matrix k(X,X) with respect to each hyperparameter
            % g = Gp.CovFn.gradient(X, GP)
            %
            % Gp.CovFn is an instantiation of the Remap covariance function class
            %
            % Inputs:   X = Input locations (D x N matrix)
            %           GP = The GP class instance can be passed to give the covariance function access to its properties
            % Outputs:  g = gradient of the covariance matrix k(X,X) with respect to each hyperparameter. Cell Array (1 x Number of Hyperparameters). Each Element is a N x N matrix


            parin = this.getCovPar(GP);
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
        function v = pointval(this, x_star, GP)
            % Efficient case for getting diag(k(X,X))
            % v = Gp.CovFn.pointval(X, GP)
            %
            % Gp.CovFn is an instantiation of the Remap covariance function class
            %
            % Inputs : X = Input locations (D x n matrix)
            % Output : v = diag(k(X,X))
        
             parin = this.getCovPar(GP);
             [D,N1] = size(x_star); %number of points in X1
             if (length(parin)~=this.npar(D))
                error('tacopig:inputInvalidLength','Wrong number of hyperparameters!');
            end
            par = parin(this.indx);
            v = this.covfn.pointval(x_star, par);
        end
        
    end
end
