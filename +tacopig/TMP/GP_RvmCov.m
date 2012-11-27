% Relevance vector machine covariance function
% This was tested using the @gaussbasisfun as the basis fun

classdef GP_RvmCov < GP_CovFunc
   
    properties
       basisfun
       numbasisfun 
    end
    
    
    methods
         function this = GP_RvmCov(basisfunptr,numbasisfun)
           this.basisfun = basisfunptr;
           this.numbasisfun = numbasisfun;
        end   
        
        function n_theta = npar(this,D)
            xsample = rand(D,1);
            numpar_basisfun = feval(this.basisfun,xsample);
            n_theta = numpar_basisfun*this.numbasisfun+this.numbasisfun+1; % params for each basis function, the weight of each basis fun and the sigma_f 
        end
        
        %Evaluate covariance between training inputs
       function K = Keval(this, X1, par)
            [D,N1] = size(X1); %number of points in X1
            xsample = rand(D,1);
            numpar_basisfun = feval(this.basisfun,xsample);
            Basis_Params = par(1:numpar_basisfun*this.numbasisfun);
            Weights = par(numpar_basisfun*this.numbasisfun+1:end-1);
            Sigma_f = par(end);
            
            K       = zeros(N1,N1);
            
            for i = 1:this.numbasisfun
              Basis = feval(this.basisfun,X1,Basis_Params(i*numpar_basisfun-(numpar_basisfun-1):i*numpar_basisfun));
              Basis_w = Basis .* Weights(i);  
              K       = K + Basis_w' * Basis;
            end
            
            K = Sigma_f*K;
            
       end
          
       %Evaluate covariance between two different sets of inputs
       function K = eval(this, X1, X2, par)
            [D,N1] = size(X1); %number of points in X1
            N2 = size(X2,2); %number of points in X2
            xsample = rand(D,1);
            numpar_basisfun = feval(this.basisfun,xsample);
            Basis_Params = par(1:numpar_basisfun*this.numbasisfun);
            Weights = par(numpar_basisfun*this.numbasisfun+1:end-1);
            Sigma_f = par(end);
            
            K       = zeros(N1,N2);
            
            for i = 1:this.numbasisfun
              BasisX1 = feval(this.basisfun,X1,Basis_Params(i*numpar_basisfun-(numpar_basisfun-1):i*numpar_basisfun));
              BasisX1_w = BasisX1 .* Weights(i); 
              BasisX2 = feval(this.basisfun,X2,Basis_Params(i*numpar_basisfun-(numpar_basisfun-1):i*numpar_basisfun));              
              K       = K + BasisX1_w' * BasisX2;
            end
            
            K = Sigma_f*K;
            
       end     
        
    end
end
