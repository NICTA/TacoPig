
%% Squared Exponential with Noisy Inputs
% Based on Girard's analytical solution
% Assumes noise covariance is a diagnal matrix
% X_var = [VarX1Dim1 VarX2Dim1 ... ; VarX1Dim2 VarX2Dim2]


classdef GP_SqExpNoisyCov < GP_CovFunc
   
    % Most covariance functions will be static
    properties
        X_var;  
    end
  
    methods
        
        function this = GP_SqExpNoisyCov(X_var)
            this.X_var = X_var;
        end
            
        function n_theta = npar(this,D)
            n_theta = D+1; % each dimension + signal variance
        end
        
        function [K] = Keval(this,X, par)
            %Assumes X1 is training and X2 is test.
            
            [D,N] = size(X); %number of points in X1

            ExpPower = zeros(N);
            IdentUncert = ones(N);
            SX1 = this.X_var;
            SX2 = SX1;
            for i = 1:D

                XX1 = repmat(X(i,:),N,1)';
                XX2 = repmat(X(i,:),N,1);
                dist = XX1-XX2;

                SSX1 = repmat(SX1(i,:),N,1)';
                SSX2 = repmat(SX2(i,:),N,1);
                AddSSX1SSX2 = SSX1+SSX2;

                UncertainLengthscale = par(i)^2+AddSSX1SSX2;

                ExpPower = ExpPower + dist.^2./UncertainLengthscale;

                IdentUncert = IdentUncert.*(1+ par(i)^-2*AddSSX1SSX2);
            end

            K = par(end)^2*IdentUncert.^(-0.5).*exp(-0.5*ExpPower);

         end
        
        
            function [K] = eval(this,X,x, par)
                %Assumes X1 is training and X2 is test.

                [D,N] = size(X); %number of points in X1
                [d,n] = size(x); %number of points in X1

            ExpPower = zeros(N,n);
            IdentUncert = ones(N,n);
            SX1 = this.X_var;
            SX2 = zeros(d,n);
            for i = 1:D

                XX1 = repmat(X(i,:),n,1)';
                XX2 = repmat(x(i,:),N,1);
                dist = XX1-XX2;

                SSX1 = repmat(SX1(i,:),n,1)';
                SSX2 = repmat(SX2(i,:),N,1);
                AddSSX1SSX2 = SSX1+SSX2;

            %     UncertainLengthscale = (1+par(i).*AddSSX1SSX2)./par(i);
                UncertainLengthscale = par(i)^2+AddSSX1SSX2;

                ExpPower = ExpPower + dist.^2./UncertainLengthscale;

                IdentUncert = IdentUncert.*(1+ par(i)^-2*AddSSX1SSX2);
            end

            K = par(end)^2*IdentUncert.^(-0.5).*exp(-0.5*ExpPower);                
                
            end
        
%         Also overload the point covariance kx*x* - its trivial
            function v = pointval(this,x_star, par)
                v = par(end).^2 * ones(1,size(x_star,2));
            end
        
    end
end
