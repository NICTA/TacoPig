classdef Exp < tacopig.covfn.CovFunc

    methods
        
        function n_theta = npar(this, D)
            n_theta = D+1; % each dimension + signal variance
        end
        
        function [K z] = eval(this, X1, X2, GP)
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
            par = this.getCovPar(GP);
            [D,N1] = size(x_star); %number of points in X1
            if (length(par)~=D+1)
                error('tacopig:inputInvalidLength','Wrong number of hyperparameters!');
            end
            v = par(end).^2 * ones(1,size(x_star,2));
        end
        
    end
end
