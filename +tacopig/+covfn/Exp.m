% The Exponential covariance function class
% k(X1,X2) = Sigma_f^2*exp(-sqrt((X1-X2)'*diag(Lengthscales.^-2)*(X1-X2)))
%
% X1 and X2 are input matrices of dimensionality D x N and D x M, respectively. 
% D is the dimesionality of the data. 
% N and M are the number of observations in the matrix X1 and X2, respectively.
% k(X1,X2) is a N x M covariance matrix
%
% Note on hyperparameters:
% The hyperparameter vector has a dimensionality of 1 x (D+1) and is of the form [Lengthscales, Sigma_f]
% Sigma_f is the signal variance hyperpameter 
% Lengthscales is a 1 x D matrix

classdef Exp < tacopig.covfn.CovFunc
      
    properties(Constant)
        ExampleUsage = 'tacopig.covfn.Exp()'; %Instance of class created for testing
    end
      
    methods

        function this = Exp() 
            % Exp covariance function constructor
            % GP.CovFn = tacopig.covfn.Mat3()
            % Gp.CovFn is an instantiation of the SqExp covariance function class 
        end    
        
        function n_theta = npar(this, D)
            % Returns the number of hyperparameters required by all of the children covariance functions
            %
            % n_theta = GP.covfunc.npar(D)
            %
            % Gp.CovFn is an instantiation of the Exp covariance function class
            % Inputs: D (Dimensionality of the dataset)
            % Outputs: n_theta (number of hyperparameters required by all of the children covariance function)
            n_theta = D+1; % each dimension + signal variance
        end
        
        function [K z] = eval(this, X1, X2, GP)
            % Get covariance matrix between input sets X1,X2 
            %
            % K = GP.covfunc.eval(X1, X2, GP)
            %
            % Gp.CovFn is an instantiation of the Exp covariance function class
            % Inputs:  X1 = D x N Input locations
            %          X2 = D x M Input locations
            %          GP = The GP class instance can be passed to give the covariance function access to its properties
            % Outputs: K = covariance matrix between input sets X1,X2 (N x M)
            
            par = this.getCovPar(GP);
            [D,N1] = size(X1); %number of points in X1
            N2 = size(X2,2); %number of points in X2
            if (D~=size(X2,1))
                error('tacopig:dimMismatch','Dimensionality of X1 and X2 must be the same.');
            end
            
            if (length(par)~=D+1)
                error('tacopig:inputInvalidLength','Wrong number of hyperparameters for NegExp');
            end
            %Compute weighted squared distances:
            w = par(1:D)'.^(-2);
            XX1 = sum(w(:,ones(1,N1)).*X1.*X1,1);
            XX2 = sum(w(:,ones(1,N2)).*X2.*X2,1);
            X1X2 = (w(:,ones(1,N1)).*X1)'*X2;
            XX1T = XX1';
            % numerical effects can drive z slightly negative 
            z = max(0,XX1T(:,ones(1,N2)) + XX2(ones(1,N1),:) - 2*X1X2);
            K = par(D+1)^2 * exp(-sqrt(z)); 
        end
        
        function [g] = gradient(this, X, GP)
            % Returns gradient of the covariance matrix k(X,X) with respect to each hyperparameter
            % g = Gp.CovFn.gradient(X, GP)
            %
            % Gp.CovFn is an instantiation of the Exp covariance function class
            %
            % Inputs:   X = Input locations (D x N matrix)
            %           GP = The GP class instance can be passed to give the covariance function access to its properties
            % Outputs:  g = gradient of the covariance matrix k(X,X) with respect to each hyperparameter. Cell Array (1 x Number of Hyperparameters). Each Element is a N x N matrix
            
            par = this.getCovPar(GP);
            %Gradient currently not working correctly.
            epsilon = 1e-10; 
            % Same as K?
            [Kg z] = this.eval(X, X, par);
            z = z+epsilon*eye(size(z)); %jitter added to avoid divide by zero
            [d,n] = size(X);
            g = cell(1,d+1);
            for i=1:d
 
            %Compute weighted squared distance
                row = X(i,:);
                XXi = row.*row;
                XTXi = row'*row;
                XXi = XXi(ones(1,n),:);
                zi = max(0,XXi+XXi'-2*XTXi); %zi is not normalised by l
                wi = par(i)^(-3);
                g{i} = par(d+1)^2.*zi.* wi.*(z).^(-0.5).*exp(-sqrt(z));
            end
            g{d+1} = Kg*(2/par(d+1));
        end
        
        
        % Also overload the point covariance kx*x* - its trivial
        function v = pointval(this, x_star, GP)
            % Efficient case for getting diag(k(X,X))
            % v = Gp.CovFn.pointval(X, GP)
            %
            % Gp.CovFn is an instantiation of the Exp covariance function class
            %
            % Inputs : X = Input locations (D x n matrix)
            % Output : v = diag(k(X,X))
            
            par = this.getCovPar(GP);
            [D,N1] = size(x_star); %number of points in X1
            if (length(par)~=D+1)
                error('tacopig:inputInvalidLength','Wrong number of hyperparameters!');
            end
            v = par(end).^2 * ones(1,size(x_star,2));
        end
        
    end
end
