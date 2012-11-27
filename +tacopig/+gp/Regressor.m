% A standard GP-regression model

classdef Regressor < tacopig.gp.GpCore
    
    properties
        mu                  % Evaluated Mean
        K                   % Evaluated Covariance Matrix
        meanpar             % Mean Function Parameters
        covpar              % Covariance Function Parameters
        noisepar            % Noise Function Parameters
        alpha               % cached K^-1y
        opts                % Optimiser options vector
        
        objective_function  % optimisation criteria for learning
        solver_function     % External optimisation toolbox function
        factorisation       % Factorisation method for inference
        
        factors             % Factors produced by the factorisation
        has_been_solved     % Flag to stop premature querying of a model
        lml                 % log marginal likelihood of training data
    end
    
   
    methods
        
        % Constructor
        function this = Regressor()
             this.opts =  optimset('LargeScale','off','MaxIter',1000,'Diagnostics', 'on',...
                 'GradObj', 'on', 'Display', 'iter','TolFun',  1e-10, 'MaxFunEvals', 5000);
             this.factorisation = 'CHOL';
             this.objective_function = @tacopig.objectivefn.NLML;
             this.solver_function = @fminunc;
             this.has_been_solved = 0;
        end
        
        % Provide convenient access to optimisation settings
        function optimset(varargin)
            this = varargin{1};
            this.opts = optimset(this.opts, varargin{2:end});
        end
        function out = optimget(varargin)
            this = varargin{1};
            out = optimget(this.opts, varargin{2:end});
        end
        
        % Run the GP Inference
        function solve(this)

            % Check validity of current configuration
            this.check(); 
            
            % Invoke the components
            mu = this.MeanFn.eval(this.X, this.meanpar);
            K0 = this.CovFn.Keval(this.X, this.covpar);
            noise = this.NoiseFn.eval(this, this.noisepar);
            
            K = K0 + noise;
            ym = (this.y - mu)';
            
            % We offer different factorisation methods
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
                this.factors.type = 'svd';
            elseif strcmpi(this.factorisation, 'chol')
                L = chol(K, 'lower');
                this.alpha = L'\(L\ym);
                this.K = K; 
                this.factors.L = L;
                this.factors.type = 'chol';
                this.mu = mu;
            else
                error('tacopig:badConfiguration','Invalid factorisation method.');    
            end
            
            % Save the solved outputs:
            this.K = K; 
            this.mu = mu;
            this.lml = -tacopig.objectivefn.NLML(this, [this.meanpar, this.covpar, this.noisepar]);
            this.has_been_solved = 1;
        end

        
        % Query the model after it has been solved
        function [mu_star, var_star, var_full] = query(this, x_star, batches)
            
            % The user can (optionally) split the data into batches)
            if (nargin<3)
                batches = 1;
            end
            if (~this.has_been_solved)
                error('tacopig:badConfiguration', 'GP must be solved first using GP.solve.');
            end
            this.check();
            
            % Get input lengths
            N = size(this.X,2); 
            nx = size(x_star,2);
            
            if abs(round(batches)-batches)>1e-16
                error('tacopig:inputInvalidType', 'Batches must be an integer');
            end
            batches = round(batches);
            % The user can also (optionally get the full variance 
            if (nargout==3)
                % provided they havent split the data
                if (batches > 1)
                    error('tacopig:badConfiguration', 'Cannot obtain full variance in batches');
                end
                
                if (nx>2000)
                    disp(['Warning: Large number of query points.'...
                        'This may result in considerable computational time.'...
                        'Press Ctrl+C to abort or any other key to continue.'])
                    pause
                end
                
            end
            
            
            mu_star = zeros(1,nx);
            if (nargout>1) 
                var_star   = zeros(1,nx);
            end
                        
            use_svd = false;
            use_chol = false;
            if strcmpi(this.factors.type, 'svd')
                factorS = sqrt(this.factors.S2);
                use_svd = true;
            elseif strcmpi(this.factors.type, 'chol')
                use_chol = true;
            else
                error('tacopig:badConfiguration', 'Factorization is not recognized.');
            end
            
            % we are currently handling the possibility of multi-task with
            % common points as a general case of GP_Std
            mu_0 = this.MeanFn.eval(x_star, this.meanpar);
            
            partitions = round(linspace(1, nx+1, batches+1));
            
            for i = 1:batches
                % Handle batches
                L = partitions(i);
                R = partitions(i+1)-1;
                LR = L:R;
                
                % Progress Bar
                if (batches>1)
                    fprintf('%d to %d...\n',L, R);
                end
                
                % Compute predictive mean
                ks = this.CovFn.eval(this.X,x_star(:,LR),this.covpar)';
                mu_star(LR) = mu_0(LR) + (ks*this.alpha)';

                if (nargout>=2)
                    % Compute predictive variance
                    var0 = this.CovFn.pointval(x_star(:,LR), this.covpar);
                    if use_svd
                        %S2 = S2(:,ones(1,size(x_star(:,LR),2)));
                        v = bsxfun(@times, factorS, (ks*this.factors.SVD_U)');
                    elseif use_chol
                        v = this.factors.L\ks';
                    else
                        error('tacopig:badConfiguration', 'Factorization not implemented');
                    end
                    var_star(LR) = max(0,var0 - sum(v.*v));
                    
                    if (nargout ==3)
                        % we also want the block
                        % Can only get here if batches is set to 1
                        
                        var0 = this.CovFn.Keval(x_star(:,LR), this.covpar);
                        var_full = var0-v'*v;
                    end
                end

            end
        end
        
        % Learn the hyperparameters        
        function theta = learn(this)
            
            % Check configuration
            this.check();

            % Pack hyperparameters into a single vector
            par0 = [this.meanpar, this.covpar, this.noisepar];
            
            % It is assumed the objective function is equivalent to fminunc
            % Otherwise errors will be thrown
            [par,fval] = this.solver_function(@(theta) this.objectfun(theta), par0, this.opts);
             
            % Unpack the parameters again to save them:   
            D = size(this.X,1);
            ncovpar = this.CovFn.npar(D);
            nmeanpar = this.MeanFn.npar(D);
            nnoisepar = this.NoiseFn.npar;
            this.meanpar = par(1:nmeanpar);
            this.covpar = par(nmeanpar+(1:ncovpar));
            this.noisepar = par(ncovpar+nmeanpar+(1:nnoisepar));
            this.lml = - tacopig.objectivefn.NLML(this,par);
        end
        
        
        % Draw samples from the prior distribution        
        function [f_star] = sampleprior(this, x_star)
            this.check();
            nx = size(x_star,2);  
            if (nx>2000)
                disp(['Warning: Large number of query points.'...
                    'This may result in considerable computational time.'...
                    'Press Ctrl+C to abort or any other key to continue.'])
                pause
            end
            Mu_0= this.MeanFn.eval(x_star, this.meanpar);
            kss = this.CovFn.eval(x_star,x_star,this.covpar);
            Lss = chol(kss+diag(diag(kss))*1e-4, 'lower');
            f_star = [Lss*randn(1,nx)']'+Mu_0;
        end
       
        % Sample from the posterior distribution        
        function [f_star] = sampleposterior(this, x_star)
            this.check();
            
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
            [mu_star var_star var_full] = this.query(x_star);
            LFull  = chol(var_full+diag(var_star)*1e-4, 'lower');
            f_star = randn(1,nx)*LFull' +mu_star;
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
                error('tacopig:badConfiguration', 'Matrix factorisation should either be SVD or CHOL\n');
            elseif ~isa(this.MeanFn,'tacopig.meanfn.MeanFunc')
                error('tacopig:badConfiguration', 'Invalid Mean Function\n');
            elseif ~isa(this.CovFn,'tacopig.covfn.CovFunc')
                error('tacopig:badConfiguration', 'Invalid Covariance Function\n');
            elseif (size(this.y,1) ~= 1)
                error('tacopig:dimMismatch', 'Y is transposed?\n');
            elseif (size(this.covpar,1) > 1)
                error('tacopig:dimMismatch', 'Covpar is transposed?\n');
            end 
            ncovpar = this.CovFn.npar(D);
            nmeanpar = this.MeanFn.npar(D);
            if (size(this.covpar,2) ~= ncovpar)
                error('tacopig:inputInvalidLength', 'covpar is the wrong length.');
            elseif (size(this.meanpar,1) > 1)
                error('tacopig:dimMismatch', 'meanpar is transposed?');
            elseif (size(this.meanpar,2) ~= nmeanpar)
                error('tacopig:inputInvalidLength', 'meanpar is the wrong length.');
            end
            return;
        end

        % Compare analytic with numerical gradients       
        function [analytic, numerical] = check_gradients(this,params)
            
            if strcmpi(optimget(this.opts, 'GradObj'), 'off')
                error('tacopig:badConfiguration', 'Gradients are turned off.');
            end
            fprintf('Checking gradients:\n');
            
            if (nargin==1)
                params = [this.meanpar, this.covpar, this.noisepar];
            end
            
            eps = 1e-6;
            [~,analytic] = this.objective_function(this,params);
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
        end
    end
end    