classdef SqExp < tacopig.covfn.CovFunc
   
    % Most covariance functions will be static
    methods(Static) 
        
        function n_theta = npar(D)
            if D <= 0
              error('tacopig:inputOutOfRange', 'Dimension cannot be < 1');
            end
            if mod(D,1) ~= 0
              error('tacopig:inputInvalidType', 'Dimension must be an integer');
            end
            n_theta = D+1; % each dimension + signal variance
        end
        
        function K = eval(X1, X2, par)
            [D,N1] = size(X1); %number of points in X1
            N2 = size(X2,2); %number of points in X2
            if D~=size(X2,1)
                error('tacopig:dimMismatch','Dimensionality of X1 and X2 must be the same');
            end
            if (length(par)~=D+1)
                error('tacopig:inputInvalidLength','Wrong number of hyperparameters for SqExp');
            end
            %Compute weighted squared distances:
            w = par(1:D)'.^(-2);
            XX1 = sum(w(:,ones(1,N1)).*X1.*X1,1);
            XX2 = sum(w(:,ones(1,N2)).*X2.*X2,1);
            X1X2 = (w(:,ones(1,N1)).*X1)'*X2;
            XX1T = XX1';
            % numerical effects can drive z slightly negative 
            z = max(0,XX1T(:,ones(1,N2)) + XX2(ones(1,N1),:) - 2*X1X2);
            K = par(D+1)^2 * exp(-0.5*z); 
        end
        
        function [g] = gradient(X, par)
            
            % Same as K?
            Kg = tacopig.covfn.SqExp.eval(X, X, par);
            
            [d,n] = size(X);
            g = cell(1,d+1);
            for i=1:d
                %Compute weighted squared distance
                row = X(i,:);
                XX = row.*row;
                XTX = row'*row;
                XX = XX(ones(1,n),:);
                z = max(0,XX+XX'-2*XTX);
                w = par(i)^(-3);
                g{i} = w*Kg.*z;
            end
            g{d+1} = Kg*(2/par(end));
        end
        
        
        % Also overload the point covariance kx*x* - its trivial
        function v = pointval(x_star, par)
            v = par(end).^2 * ones(1,size(x_star,2));
        end
        
    end
end
