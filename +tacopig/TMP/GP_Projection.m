% Allows a projection matrix to relate the observations to the underlying kernel
classdef GP_Projection < GP_CovFunc
    
    properties
       G
       covfn
       %maxblock
    end
    
    methods
        function this = GP_Projection(G, covfn)
           if ~isa(covfn,'GP_CovFunc')
               error([class(this),': must specify a valid covariance function']); 
           end
           this.G = G;
           this.covfn = covfn;
        end   
            
        function n_theta = npar(this,D)
            n_theta = this.covfn.npar(D);
        end

        function K = Keval(this, X, par)
            eps = 1e-4;
            Gi = this.G;
            K0 = this.covfn.Keval(X, par);
            K0 = K0 + eps*min(diag(K0))*eye(size(K0));
            K = Gi*K0*Gi';
            K = K + eps*min(diag(K))*eye(size(K));
        end
        
        function g = gradient(this,X, par)
            % Also fix this up to use large blocks as neccessary...
            GI = this.G;
            g0 = this.covfn.gradient(X, par);
            g = g0;
            for i=1:length(par)
                g{i} = GI*g{i}*GI';
            end
        end
                
        function v = pointval(this, x_star, par)
            % G is not involved within the projection
            v = this.covfn.pointval(x_star, par);
        end
        
        function K = eval(this, X1, X2, par)
            K = (this.G)*this.covfn.eval(X1, X2, par);
        end
    end
end

           % if (true) %(nx<=this.maxblock) % split into blocks of K... makes learning very costly
            
%             else
%                 % What size will K be
%                 nk = size(G,1);
%                 K = zeros(nk);
%                 step = this.maxblock;
%                 for i=1:step:nx
%                     
%                     
%                     LR = i:min(i+step-1,nx);
%                     Gi = G(:, LR);
%                     
%                     K0 = this.covfn.Keval(X(:,LR), par);
%                    
% %                     size(Gi)
% %                     size(K0)
%                     K = K + Gi*K0*Gi';
%                     
%                     for j=(i+step):step:nx %
% 
%                         TB = j:min(j+step-1,nx);
%                         Gj = G(:, TB);
%                         
%                         K0 = this.covfn.eval(X(:,TB),X(:,LR),par);
% 
%                         % can be nonsquare but must be square???
%                         K = K + Gi*K0'*Gj' + Gj*K0*Gi'; % The Kij case and Kji case
%                     end
%                 end
%             end

