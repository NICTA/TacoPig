classdef GP_SqExpNonStatCov < GP_CovFunc

    % The non-stationary covariance function described in Plagemann ECML08
    
    properties
       lengthfunptr
       lengthfuncomplexity % a vector passed to the lengthscale function to describe its complex e.g. the degree in the case of a polynomial function
    end
    
    methods
        
        function this = GP_SqExpNonStatCov(lengthfunptr,lengthfuncomplexity)      
           this.lengthfunptr = lengthfunptr; %pointer to the lengthscale function
           this.lengthfuncomplexity = lengthfuncomplexity; %pointer to the lengthscale function
        end  
        
        function n_theta = npar(this,D)
            xsample = rand(D,1);
            numpar_lengthfun = feval(this.lengthfunptr,xsample,this.lengthfuncomplexity);
            n_theta = numpar_lengthfun+1; % each dimension + signal variance
        end
        
        function K = eval(this, X1, X2, par)
            [D,N1] = size(X1); %number of points in X1
            N2 = size(X2,2); %number of points in X2
            
            LengthFunPar = par(1:end-1);
            Sigma1 = abs(feval(this.lengthfunptr,X1,this.lengthfuncomplexity,LengthFunPar));
            Sigma2 = abs(feval(this.lengthfunptr,X2,this.lengthfuncomplexity,LengthFunPar));
            Sig1SquaredRep = (Sigma1(ones(1,N2),:)').^2;
            Sig2SquaredRep = (Sigma2(ones(1,N1),:)).^2;
            SigMat = 0.5*Sig1SquaredRep+0.5*Sig2SquaredRep;
%             [xtrain ytrain]=meshgrid(-30:4:30,-30:4:30);
%             figure; surf(xtrain,ytrain,reshape(Sigma1,size(xtrain)))
            
            X1Stacked = reshape(X1',size(X1,2),1,size(X1,1));
            X1rep = X1Stacked(:,ones(1,N2),:);

            X2Stacked = reshape(X2',1, size(X2,2),size(X2,1));
            X2rep = X2Stacked(ones(1,N1),:,:);

            Dist = X1rep-X2rep;
            K = ones(N1,N2);
            for i=1:D
                K =K.*(Sig1SquaredRep.^0.25).*(Sig2SquaredRep.^0.25)...
                .*(SigMat).^(-0.5).*exp(-Dist(:,:,i).^2./SigMat);
            end
            K = par(end)^2*K;
        end
  
        
    end
end
