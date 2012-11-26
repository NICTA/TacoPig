% A standard GP-regression model

classdef GP_STD < GP_Model
    
    properties
        mu              % Evaluated Mean
        K               % Evaluated Covariance Matrix
        meanpar         % Mean Function Parameters
        covpar          % Covariance Function Parameters
        noisepar        % Noise Function Parameters
        alpha           % cached K^-1y
        opts            % Optimiser options vector
        
        factorisation   % Factorisation method for inference
        factors         % Factors produced by the factorisation
        has_been_solved 
        
        lml
        objective_function
        solver_function
        
    end
    
    methods
        
        % Constructor
        function this = GP_STD()
             this.opts =  optimset('LargeScale','off','MaxIter',1000,'Diagnostics', 'on',...
                 'GradObj', 'on', 'Display', 'iter','TolFun',  1e-10, 'MaxFunEvals', 5000);
             this.factorisation = 'SVD';
             this.objective_function = @GP_LMLG_FN;
             this.solver_function = @fminunc;
             this.has_been_solved = 0;
        end
        
% Access to optimisation settings
        function optimset(varargin)
            this = varargin{1};
            this.opts = optimset(this.opts, varargin{2:end});
        end
        function out = optimget(varargin)
            this = varargin{1};
            out = optimget(this.opts, varargin{2:end});
        end
        
        
        
% Solve for K, lml, etc.          
        function solve(this)
            this.check();
            mu = this.MeanFn.eval_y(this.X, this.meanpar);
            K0 = this.CovFn.Keval(this.X, this.covpar);
            noise = this.NoiseFn.eval(this, this.noisepar);
            ym = (this.y - mu)';
            K = K0 + noise;
            
            if strcmpi(this.factorisation, 'svd')
                [U,S,V] = svd(K);
                S2 = diag(S);
                S2(S2>0) = 1./S2(S2>0);
                invK = V*diag(S2)*U';
                this.alpha = invK*ym;
                this.factors.S2 = S2;
                this.factors.SVD_U = U;
                this.factors.SVD_S = S;
                this.factors.SVD_V = V;
            elseif strcmpi(this.factorisation, 'chol')
                L = chol(K, 'lower');
                this.alpha = L'\(L\ym);
                this.K = K; 
                this.factors.L = L;
                this.mu = mu;
            else
                error('Invalid factorisation!');    
            end
            
            this.K = K; 
            this.mu = mu;
            this.lml = -GP_LMLG_FN(this, [this.meanpar, this.covpar, this.noisepar]);
            this.has_been_solved = 1;
        end
        
        
% Query the model               
        function [mu_star, var_star, var_full] = query(this, x_star)
            if (~this.has_been_solved)
                error('GP must be solved first using GP.solve.');
            end
            
            this.check();
                            
            N = size(this.X,2); 
            nx = size(x_star,2);

            array_limit = 10*1024*1024/8; % number of doubles in 10mb
            
            if (N*nx>array_limit)
                large_query = true;
                fprintf('Evaluating query points in 10mb batches:\n')
            else
                large_query = false;
            end
            
            mu_star   = zeros(1,nx);
            if (nargout>1) 
                var_star   = zeros(1,nx);
            end
            ee   = 0; % loop counter

            step = floor(array_limit/N);
            if (step==0)
                warning('X is very large - querying single file.\n');
                step = 1;
            end
                        
            use_svd = false;
            use_chol = false;
            if strcmpi(this.factorisation, 'svd')
                % Cache - common to all queries
                S2 = sqrt(this.factors.S2);
                use_svd = true;
            elseif strcmpi(this.factorisation, 'chol')
                use_chol = true;
            else
                error('Invalid choice of factorisation method');
            end
            
            % we are currently handling the possibility of multi-task with
            % common points as a general case of GP_Std
            mu_0 = this.MeanFn.eval(x_star, this.meanpar);
            uptomu0 = 1;
            for i = 1:step:nx
                % progress bar
                if (large_query)
                    fprintf('%d...',i);
                end
                LR = i:min(nx, i+step-1);
                ks = this.CovFn.eval(this.X,x_star(:,LR),this.covpar)';
                product = (ks*this.alpha)';
                
                nlr = size(LR,2);
                nrep = size(ks,1)/nlr; % repetition of the input data...
                mu_star(1:nrep,LR) = reshape(product, [nlr,nrep])';
                mu_1 = reshape(mu_0, size(mu_star'))';
                for j=1:nrep
                    mu_star(j,LR) = mu_star(j,LR) + mu_1(j,LR);
                end
                % Realised a bug that mis-calculated the variance in this
                % case...
                if (nargout==2) 
                    var0 = this.CovFn.pointval(x_star(:,LR), this.covpar);
                    if use_svd
                        %S2 = S2(:,ones(1,size(x_star(:,LR),2)));
                        v = bsxfun(@times, S2, (ks*this.factors.SVD_U)');
                    elseif use_chol
                        v = this.factors.L\ks';
                    end
                    predvar = max(0,var0 - sum(v.*v));
                    var_star(1:nrep,LR) = reshape(predvar, [nlr,nrep])'; % is that more correct?
                end

            end
            fprintf('\n');
            if (nargout ==3)
                     var0 = this.CovFn.eval(x_star, x_star, this.covpar);
                    if use_svd
                        v = bsxfun(@times, S2, (ks*this.factors.SVD_U)');
                    elseif use_chol
                        v = this.factors.L\ks';
                    end
                    var_star = max(0,diag(var0)' - sum(v.*v));
                    var_full = var0-v'*v;
            end
        end
          
        
% Learn the hyperparameters        
        function theta = learn(this)
            if (~this.check())
                error('Abort.');
            end
            
            %Matlab optimisation toolbox
            par0 = [this.meanpar, this.covpar, this.noisepar];
            [par,fval] = this.solver_function(@(theta) this.objectfun(theta), par0, this.opts);
             
            % Unpack the parameters again to save them!   
            D = size(this.X,1);
            ncovpar = this.CovFn.npar(D);
            nmeanpar = this.MeanFn.npar(D);
            nnoisepar = this.NoiseFn.npar;
            this.meanpar = par(1:nmeanpar);
            this.covpar = par(nmeanpar+(1:ncovpar));
            this.noisepar = par(ncovpar+nmeanpar+(1:nnoisepar));
            this.lml = -GP_LMLG_FN(this,par);
        end
        
% Sample from the prior distribution        
        function [f_star] = sampleprior(this, x_star)
            if (~this.check())
                error('Abort.');
            end
            nx = size(x_star,2);  
            if nx>20000
                disp(['Warning: Large number of query points.'...
                    'This may result in considerable computational time.'...
                    'Press Ctrl+C to abort or any other key to continue.'])
                pause
            end
            f_star   = zeros(1,nx);
            Mu_0= this.MeanFn.eval(x_star, this.meanpar);
            kss = this.CovFn.eval(x_star,x_star,this.covpar);
            Lss = chol(kss+eye(size(kss))*1e-10, 'lower');
            f_star = [Lss*randn(1,nx)']'+Mu_0;
        end
       
% Sample from the posterior distribution        
        function [f_star] = sampleposterior(this, x_star)
            if (~this.check())
                error('Abort.');
            end
            
            if (~this.has_been_solved)
                error('GP must be solved first using GP.solve().');
            end
            
            nx = size(x_star,2);  
            N  = size(this.X,2);
            if (nx+N)>20000
                disp(['Warning: Large number of query points.'...
                    'This may result in considerable computational time.'...
                    'Press Ctrl+C to abort or any other key to continue.'])
                pause
            end
            f_star = zeros(1,nx);
            [mu_star var_star var_full] = this.query(x_star);
            LFull = chol(var_full+eye(size(var_full))*1e-8, 'lower');
            fFull = [LFull*randn(1,nx)']';
            f_star = fFull+mu_star;

        end
        
        % Pass through objective function
        function [objective, objective_grad] = objectfun(this, parvec)
            if (nargout == 2)
                [objective, objective_grad] = this.objective_function(this, parvec);
            elseif (nargout == 1)
                 objective = this.objective_function(this, parvec);
            else
                error('Wrong number of output arguments');
            end
        end
        
        
        % Check the inputs are of an appropriate form        
        function check(this)
            % Check that all the values are properly set
            [D,N] = size(this.X);
            use_svd = strcmpi(this.factorisation, 'svd');
            use_chol = strcmpi(this.factorisation, 'chol');
            if ((use_svd==0)&&(use_chol==0))
                error('Matrix factorisation should either be SVD or CHOL\n');
            elseif ~isa(this.MeanFn,'GP_MeanFunc')
                error('Not a valid Mean Function\n');
            elseif ~isa(this.CovFn,'GP_CovFunc')
                error('Not a valid Covariance Function\n');
            elseif (size(this.y,1) ~= 1)
                error('Y is transposed?\n');
            elseif (size(this.covpar,1) > 1)
                error('covpar is transposed?\n');
            end 
            ncovpar = this.CovFn.npar(D);
            nmeanpar = this.MeanFn.npar(D);
            if (size(this.covpar,2) ~= ncovpar)
                error('covpar is the wrong length.');
            elseif (size(this.meanpar,1) > 1)
                error('meanpar is transposed?');
            elseif (size(this.meanpar,2) ~= nmeanpar)
                error('meanpar is the wrong length.');
            end
            return;
        end

        % Compare analytic with numerical gradients       
        function [analytic, numerical] = check_gradients(this,params)
            
            if strcmpi(optimget(this.opts, 'GradObj'), 'off')
                fprintf('You have gradients turned off...\n');
                return;
            end
            fprintf('Checking gradients:\n');
            
        
            if (nargin==1)
                params = [this.meanpar, this.covpar, this.noisepar];
            end
            
            eps = 1e-6;
            [val0,analytic] = this.objective_function(this,params);
            numerical = zeros(size(params));
            for i = 1:length(params)
               params2 = params;
               params3 = params;
               params2(i) = params2(i) + eps;
               params3(i) = params3(i) - eps;
               val1 = this.objective_function (this,params2);
               val2 = this.objective_function (this,params3);
               numerical(i) = (val1-val2)/(2*eps);
            end
            analytic = analytic';
            fprintf('Analytic (top), Numerical (bottom)\n');
            disp([analytic;numerical]);
            
            %disp([analytic; numerical]);
            
        end
    end
end    