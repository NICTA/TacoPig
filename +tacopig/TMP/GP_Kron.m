% Repeats a single (weighted) covariance matrix over multiple tasks
% This is a straightforward (but fairly inflexible) way to approach
% multi-task modeling. It also requires the tasks to be sampled on the same
% X. The general case of different X for each task is coming next....
% Alistair 8/Aug/2012

classdef GP_Kron < GP_CovFunc

    properties
       ntask
       covfn
    end

    methods
        function this = GP_Kron(nt, covfn)
           if ~isa(covfn,'GP_CovFunc')
               error([class(this),': must specify a valid covariance function']); 
           end
           this.ntask = nt;
           this.covfn = covfn;
        end   
            
        function n_theta = npar(this,D)
            nt = this.ntask;
            n_theta = this.covfn.npar(D) + nt*(nt+1)/2;
        end

        
        % This function unpacks the hyperparameters into a PSD signal
        % variance matrix, and a vector for the underlying function
        function [W, pcov] = split(this,par)
            nt = this.ntask;
            taskh = nt*(nt+1)/2;
            ptask = par(1:taskh);
            pcov = par(taskh+1:end);
            W = tril(ones(nt));
            W(W>0) = ptask;
            W = W*W'; % LU
        end
        
        function K = Keval(this, X, par)
            [W, pcov] = this.split(par);
            K = kron(W, this.covfn.Keval(X, pcov));
        end

        function g = gradient(this,X, par)
             nt = this.ntask;
             taskh = nt*(nt+1)/2;
             ptask = par(1:taskh);
             pcov = par(taskh+1:end);
             
             W = tril(ones(nt));
             W(W>0) = ptask;
             W = W*W'; % LU
             
             K0 = this.covfn.Keval(X, pcov);
             G0 = this.covfn.gradient(X, pcov);
             %G0 = G0(1:end-1);
             %K = kron(W, K0);
             
             G1 = cell(1,taskh);
             
             % The first taskh hypers control the sigvar W
             % i and j coordinates of where the hypers are inserted into L
             L = tril(ones(nt));
             [jj,ii] = meshgrid(1:nt, 1:nt);
             sampi = ii(L>0);
             sampj = jj(L>0);
             L(L>0) = ptask;
             % Now the derivative of W=LL' wrt the i,j th element of L is
             % equal to the jth column of L added to the ith row and the
             % ith column of zeros(nt) !!!
             for s = 1:taskh
                 i = sampi(s);
                 j = sampj(s);
                 col = L(:,j);
                 deriv = zeros(nt);
                 deriv(:,i) = col;
                 deriv(i,:) = col;
                 deriv(i,i) = deriv(i,i)*2;
                 %derivs{s} = deriv;
                 
                 G1{s} = kron(deriv, K0);
             end
             
             % And then the much simpler task of getting the gradient wrt
             % the underlying function
             
             for i=1:length(G0)
                 G0{i} = kron(W, G0{i});
             end
             
             % Then join the two together: 
             g = {G1{:},G0{:}};
         end
                
        % x_star won't be large thanks to batch calling from the GP
        % So we don't need to be as worried about pointval and eval
        function v = pointval(this, x_star, par)
            %v = this.covfn.pointval(x_star, par);
            [W, pcov] = this.split(par);
            v = kron(diag(W)', this.covfn.pointval(x_star, pcov));
        end
        
        function K = eval(this, X1, X2, par)
            [W, pcov] = this.split(par);
            K = kron(W, this.covfn.eval(X1, X2, pcov));
        end
        
    end
end
