 function [crossvalerr crossvalerr_grad] = GP_CRSVAL_FN(this, parvec)
            
            % Get configuration
            use_svd = strcmpi(this.factorisation, 'svd');
            use_chol = strcmpi(this.factorisation, 'chol');
            

            % Unpack the parameters
            [D,N] = size(this.X);
            ncovpar = this.CovFn.npar(D);
            nmeanpar = this.MeanFn.npar(D);
            nnoisepar = this.NoiseFn.npar;
            meanpar = parvec(1:nmeanpar);
            covpar = parvec(nmeanpar+1:nmeanpar+ncovpar);
            noisepar = parvec(nmeanpar+ncovpar+1:nmeanpar+ncovpar+nnoisepar);

            % save so we can stop and continue optimisation
            this.meanpar = meanpar;
            this.covpar = meanpar;
            this.noisepar = noisepar;
%             this.tmphypers = tmph;
            
            X0 = this.X;
            y0 = this.y;
            
            %Partition the data into folds
            NumFolds= 5; %Number of folds
            rand('seed',100); %Ensures same folds used all the time
            
            % Indices = crossvalind('Kfold', N, NumFolds)'; %Assigns data to folds
            % A: Turns out crossvalind is part of a toolbox - using
            % randperm instead
            Indices = mod(randperm(N), NumFolds)'+1;
            
            for i = 1:max(Indices)
                % Use all but 1 fold for training
                this.X = X0(:,Indices~=i);
                this.y = y0(:,Indices~=i);
                
                % Use remaining fold to query
                x_star = X0(:,Indices==i);
                y_star = y0(:,Indices==i);
     
                % Get mu, K and (y-mu) and gradients if needed
                mu = this.MeanFn.eval(this.X, meanpar); 
            
                K = this.CovFn.Keval(this.X, covpar);
                noise = this.NoiseFn.eval(this, noisepar);

                K = K + noise;

                ym = this.y-mu;
            
                mu_0 = this.MeanFn.eval(x_star, meanpar);
                ks = this.CovFn.eval(this.X,x_star,covpar)';
            
                % Factorise K
                if (use_svd)
                    [U,S,V] = svd(K);
                    S0 = diag(S);
                    S2 = S0;
                    S2(S2>0) = 1./S2(S2>0);
                    invK = V*diag(S2)*U';
                    alpha = (invK*ym');
                elseif (use_chol)
                    L = chol(K, 'lower');
                    alpha = L'\(L\ym');
                else
                    error('Invalid factorisation choice');
                end

                %Calc predictive mean
                mu_star = mu_0 + (ks*alpha)';
               
                %Calc predictive variance 
                var0 = this.CovFn.pointval(x_star, covpar);
                if use_svd
                    v = bsxfun(@times, sqrt(S2), (ks*U)');
                elseif use_chol
                    v = L\ks';
                end
                var_star = var0 - sum(v.*v);  
                
                % Calc neg log prob
                nlogprob(i) = sum(0.5*log(var_star)+((y_star-mu_star).^2)./(2*var_star)+0.5*log(2*pi));   
            end
            this.X = X0;
            this.y =y0;
            
            crossvalerr = sum(nlogprob); %return cross val error
            
            
            %Numerically Compute the gradient for covariance hyperparameters
            if nargout>1
                
                 val0 = GP_CRSVAL_FN (this, parvec);
                 eps = 1e-9;
                 crossvalerr_grad = zeros(size(parvec))';
                 
                for i = 1:length(parvec)
                   params2 = parvec;
                   params2(i) = params2(i) + eps;
                   val1 = GP_CRSVAL_FN (this, params2);
                   crossvalerr_grad(i) = (val1-val0)/eps;
                end
            end

        end
        