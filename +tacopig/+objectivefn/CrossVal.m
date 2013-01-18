% The cross validation objective function.
% It withholds a prespecified portion of the data (default is 1/5) and
% trains on the remainder. The performance of the model with the parameters
% is evaluated on the withheld set. This is repeated with a different portion of the data withheld 
% each time until each data point has been withheld once. The negative log likelihood
% from each trial is summed and passed as the output - crossvalerr.
%
% [nlml, nlmlg] = GP.objectivefn.NLML(parvec,NumFolds)
%
% Inputs: parvec = a concatinated vector of the mean function, covariance function and noise funtion's parameters/hyperparameters.
%         NumFolds = The number of cross validation folds
% Outputs: crossvalerr = cross validation error
%          crossvalerr_grad = the gradient of the cross validation error with respect to each of the parameters in parvec.
 
function [crossvalerr crossvalerr_grad] = CrossVal(this, parvec, NumFolds)
            
            if (nargin<3)
                NumFolds = 5;
            end
 
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
            this.covpar = covpar;
            this.noisepar = noisepar;
%             this.tmphypers = tmph;
            
            X0 = this.X;
            y0 = this.y;
            
            %Partition the data into folds
            
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
                mu = this.MeanFn.eval(this.X, this); 
                K = this.CovFn.Keval(this.X, this);
                noise = this.NoiseFn.eval(this.X, this);

                K = K + noise;
                ym = this.y-mu;
            
                mu_0 = this.MeanFn.eval(x_star, this);
                ks = this.CovFn.eval(this.X, x_star, this)';
            
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
                    error('tacopig:inputOutOfRange','Invalid factorisation choice');
                end

                %Calc predictive mean
                mu_star = mu_0 + (ks*alpha)';
               
                %Calc predictive variance 
                var0 = this.CovFn.pointval(x_star, this);
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
                
                 %val0 = tacopig.objectivefn.CrossVal(this, parvec);
                 val0 = crossvalerr;
                 eps = 1e-9;
                 crossvalerr_grad = zeros(length(parvec),1);
                 
                for i = 1:length(parvec)
                   params2 = parvec;
                   params2(i) = params2(i) + eps;
                   val1 = tacopig.objectivefn.CrossVal(this, params2);
                   crossvalerr_grad(i) = (val1-val0)/eps;
                end
            end

        end
        