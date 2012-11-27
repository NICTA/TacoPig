classdef GP_PSparseCov < GP_CovFunc
    
    % Most covariance functions will be static
    methods(Static) 
        
        function n_theta = npar(D)
            n_theta = D+1; % each dimension + signal variance
        end
        
        function [K z] = eval(X1, X2, par)
            [D,N1] = size(X1); %number of points in X1
            N2 = size(X2,2); %number of points in X2
            if D~=size(X2,1)
                error('Dimensionality of X1 and X2 must be the same');
            end
            if (length(par)~=D+1)
                error('Wrong number of hyperparameters for GP_PSparseCov');
            end
            
            %Compute weighted squared distances:
            w = par(1:D)'.^(-2);
            XX1 = sum(w(:,ones(1,N1)).*X1.*X1,1);
            XX2 = sum(w(:,ones(1,N2)).*X2.*X2,1);
            X1X2 = (w(:,ones(1,N1)).*X1)'*X2;
            XX1T = XX1';
            
            % numerical effects can drive z slightly negative 
            z = max(0,XX1T(:,ones(1,N2)) + XX2(ones(1,N1),:) - 2*X1X2);
            
            % turn into a single distance z
            z = sqrt(z)/sqrt(3);
            
            % Evaluate for a 3d problem
            j = 2.5; 
            
            z = min(1,z); % avoid imaginary style complexity
            K = par(end).^2 *max(0, (1-z).^(j+1).*((j+1)*z + 1));
            
        end
               
        % havent implemented gradients yet? 
        
        % Also overload the point covariance kx*x* - its trivial
        function v = pointval(x_star, par)
            v = par(end).^2 * ones(1,size(x_star,2));
        end
        
    end
end
