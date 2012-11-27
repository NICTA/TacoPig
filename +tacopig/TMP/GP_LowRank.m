% A low-rank Gaussian Process model that computes over a subset of regressor
% inputs (that need not belong to the training data)
% Enables inference over large datasets
% Behaviour away from training data may not be consistent with prior.

classdef GP_LowRank < GP_Model
    
    properties
        
        % Inherits the basic MeanFn, CovFn, NoiseFn, X, y
        
        % We introduce additional terms:
        XI      % induced X
        KI      % induced K
        KIX     % Cached term
        Km      % Pseudo-covariance
        palpha  % pseudo-alpha over induced points
        
        gcount
        mu
        meanpar
        covpar
        noisepar
        alpha
        spit_out_gradients
        opts
        factorisation
        factors
        invK
        lml
        objective_function
        solver_function
        has_been_solved
    end
    
    methods
        
% Construction
        function this = GP_SR()
             this.opts =  optimset('LargeScale','off','MaxIter',1000,'Diagnostics', 'on',...
                 'GradObj', 'on', 'Display', 'iter','TolFun',  1e-10, 'MaxFunEvals', 5000);
             this.factorisation = 'SVD';
             this.objective_function = @GP_SR_LMLG;
             this.solver_function = @fminunc;
             this.has_been_solved = 0;
             this.gcount = 0;
             this.spit_out_gradients = false;
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
            if (~this.check())
                error('Abort.');
            end
            
            mu = this.MeanFn.eval_y(this.X, this.meanpar);
            KI  = this.CovFn.Keval(this.XI, this.covpar);
            KIX = this.CovFn.eval(this.XI, this.X, this.covpar);
            
            % add a tiny noise to KI to keep positive definiteness
            eps = 1e-6*sum(diag(KI)); % or could use min etc
            KI  = KI + eps*eye(size(KI));
            
            % is noise a scalar or a function?
            noise = this.noisepar;% this.NoiseFn.eval(this, this.noisepar);
            Km  = (noise^2)*KI + KIX*KIX';
            ym = (this.y - mu)';
            
            %K = K + noise;
            % Now we invert Km 
            pseudoY = (KIX*ym);
            if strcmpi(this.factorisation, 'svd')
                [U,S,V] = svd(Km);
                S2 = diag(S);
                S2(S2>0) = 1./S2(S2>0);
                invK = V*diag(S2)*U';
                this.palpha = (invK*pseudoY); 
                this.invK = invK;
                this.factors.S2 = S2;
                this.factors.SVD_U = U;
                this.factors.SVD_S = S;
                this.factors.SVD_V = V;
            elseif strcmpi(this.factorisation, 'chol')
                L = chol(Km, 'lower');
                % this.palpha = L'\(L\pseudoY);
                this.invK = L'\(L\eye(size(Km)));
                this.palpha = (this.invK*pseudoY); 
                this.factors.L = L;
            else
                error('Invalid factorisation!');    
            end
            this.Km = Km; 
            this.mu = mu;
            this.lml = -GP_SR_LMLG(this, [this.meanpar, this.covpar, this.noisepar]);
            this.has_been_solved = 1;
        end
        
        
% Query the model               
        function [mu_star, var_star, var_full] = query(this, x_star)
            if (~this.has_been_solved)
                error('GP must be solved first using GP.solve.');
            end
            if (~this.check())
                error('Abort.');
            end
            
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
                ks = this.CovFn.eval(this.XI,x_star(:,LR),this.covpar)';
                product = (ks*this.palpha)';
                
                
                
                nlr = size(LR,2);
                nrep = size(ks,1)/nlr; % repetition of the input data...
                mu_star(1:nrep,LR) = reshape(product, [nlr,nrep])';
                mu_1 = reshape(mu_0, size(mu_star'))';
                for j=1:nrep
                    mu_star(j,LR) = mu_star(j,LR) + mu_1(j,LR);
                end
                if (nargout==2) 
                    vs = (this.noisepar^2)*sum((ks.*(ks*this.invK))',1);
                    predvar = max(0, vs + this.noisepar^2);%var0 - sum(v.*v));
                    var_star(1:nrep,LR) = reshape(predvar, [nrep,nlr]);
                end
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
            nnoisepar = 1;%this.NoiseFn.npar;
            this.meanpar = par(1:nmeanpar);
            this.covpar = par(nmeanpar+(1:ncovpar));
            this.noisepar = par(ncovpar+nmeanpar+(1:nnoisepar));
            this.lml = -GP_SR_LMLG(this,par);
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
            f_star = zeros(1, nx);
            [mu_star var_star var_full] = this.query(x_star);
            LFull = chol(var_full+eye(size(var_full))*1e-8, 'lower');
            fFull = [LFull*randn(1, nx)']';
            f_star = fFull+mu_star;

        end
        
        
        
        function [objective, objectiveg] = objectfun(this, parvec)
            if (nargout == 2)
                [objective, objectiveg] = this.objective_function(this, parvec);
            elseif (nargout == 1)
                 objective = this.objective_function(this, parvec);
            else
                error('Wrong number of output arguments');
            end
        end
        
        
% Check the inputs are of an appropriate form        
        function valid = check(this)
            % Check that all the values are properly set
            [D,N] = size(this.X);
            valid = false;
            
            use_svd = strcmpi(this.factorisation, 'svd');
            use_chol = strcmpi(this.factorisation, 'chol');
            if ((use_svd==0)&&(use_chol==0))
                fprintf('Matrix factorisation should either be SVD or CHOL\n');
                drawnow; 
                return;
            end
            if ~isa(this.MeanFn,'GP_MeanFunc')
                fprintf('Not a valid Mean Function\n');
                drawnow; 
                return;
            end
            if ~isa(this.CovFn,'GP_CovFunc')
                fprintf('Not a valid Covariance Function\n');
                drawnow; 
                return;
            end
            if (size(this.y,1) ~= 1)
                fprintf('Y is transposed?\n');
                drawnow; 
                return
            end
            if (size(this.covpar,1) > 1)
                fprintf('covpar is transposed?\n');
                drawnow; 
                return
            end 
            ncovpar = this.CovFn.npar(D);
            if (size(this.covpar,2) ~= ncovpar)
                error('covpar is the wrong length.');
            end
            if (size(this.meanpar,1) > 1)
                error('meanpar is transposed?');
            end 
            nmeanpar = this.MeanFn.npar(D);
            if (size(this.meanpar,2) ~= nmeanpar)
                error('meanpar is the wrong length.');
            end
            
            % No problems detected
            valid = true;
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