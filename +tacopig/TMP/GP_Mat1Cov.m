classdef GP_Mat1Cov < GP_CovFunc
   
    % Most covariance functions will be static
    methods(Static) 
        
        function n_theta = npar(D)
            n_theta = D+1; % each dimension + signal variance
        end
        
        function [K z] = eval(X1, X2, par)
            N1 = size(X1,2); %number of points in X1
            N2 = size(X2,2); %number of points in X2

            if size(X1,1)~=size(X2,1)
                error('Dimensionality of X1 and X2 must be the same');
            end
            if iscell(par)
                par = cell2num(par);
            end

            %Compute the weighted squared distances between X1 and X2
            %in a vectorised way 
            w = par(1:end-1)'.^(-2);
            XX1 = sum(w(:,ones(1,N1)).*X1.*X1,1);
            XX2 = sum(w(:,ones(1,N2)).*X2.*X2,1);
            X1X2 = (w(:,ones(1,N1)).*X1)'*X2;
            XX1T = XX1';
            z = XX1T(:,ones(1,N2)) + XX2(ones(1,N1),:) - 2*X1X2;

            K = exp(log((par(end))^2)-0.5*sqrt(z)); %put the sf2 inside the exponential
            %to avoid numerical problems.
        end
        
        function [g] = gradient(X, par)
            
          %Gradient currently not working correctly.
            % Same as K?
            % Same as K?
            Kg = GP_Mat1Cov.eval(X, X, par);
            
            [d,n] = size(X);
            g = cell(1,d+1);
            for i=1:d
                %Compute weighted squared distance
                row = X(i,:);
                XX = row.*row;
                XTX = row'*row;
                XX = XX(ones(1,n),:);
                z = max(0,XX+XX'-2*XTX);
                w = par(i)^(-2);
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
