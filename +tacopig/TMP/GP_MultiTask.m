% This is the 

classdef GP_MultiTask < GP_Model
    
    % Make properties accessible 
    properties
        % Hyperparameters
        meanpar
        covpar
        sigpar
        noisepar
        
        % Prior Parameters
        mu
        K
        SF
        alpha
        
        % Cached Posterior Parameters
        factors
        invK
        lml
        
        
        % Configuration
        ntask
        spit_out_gradients
        opts
        factorisation
        objective_function
        solver_function
        has_been_solved
    end

    
    methods

    % Construction
        function this = GP_MultiTask(ntask)
             this.opts =  optimset('LargeScale','off','MaxIter',1000,'Diagnostics', 'on',...
                 'GradObj', 'on', 'Display', 'iter','TolFun',  1e-10, 'MaxFunEvals', 5000);
             this.factorisation = 'CHOL';
             this.ntask = ntask;
             this.objective_function = @GP_MT_LMLG_FN;
             this.solver_function = @fminunc;
             this.has_been_solved = 0;
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



    % Solve the joint model in preparation for querying
        function solve(this)

            if (~this.check())
                error('Abort.');
            end

            nt = this.ntask;
            
            % Get a signal variance matrix for this case:
            SF = zeros(nt);
            SF(tril(ones(nt))>0.5) = this.sigpar(:);
            SF = SF*SF';

            X = this.X;

            % build taskmap and nx
            nx = 0;
            for t = 1:nt
                nx = nx + size(X{t},2);
            end
            taskmap = zeros(1,nx);
            upto = 0;
            for t = 1:nt
                nthis = size(X{t},2);
                taskmap(upto+(1:nthis)) = t;
                upto = upto + nthis;
            end
            K = zeros(nx); % Allocate Space
            mu = zeros(1,nx);
            covpars = this.covpar;

            for i=1:nt
                thistask = (taskmap==i);
                xi = X{i};
                mu(thistask) = this.MeanFn{i}.eval_y(xi, this.meanpar{i});
                pari = covpars{i};
                K0 = this.CovFn.eval(xi, xi, pari, pari);
                noise = this.NoiseFn{i}.eval(xi, this.noisepar{i});
                K(thistask,thistask) = SF(i,i)*K0 + noise;
                for j = i+1:nt  % dont overwrite the custom noisy one!
                    thattask = (taskmap==j);
                    xj = X{j};
                    parj = covpars{j};
                    K0 = SF(i,j)*this.CovFn.eval(xi, xj, pari, parj);
                    K(thistask,thattask) = K0;
                    K(thattask,thistask) = K0';
                end
            end   
            ym = (cat(2,this.y{:}) - mu)';

            if strcmpi(this.factorisation, 'svd')
                [U,S,V] = svd(K);
                S2 = diag(S);
                S2(S2>0) = 1./S2(S2>0);
                invK = V*diag(S2)*U';
                this.alpha = invK*ym;
                this.invK = invK;
                this.factors.S2 = S2;
                this.factors.SVD_U = U;
                this.factors.SVD_S = S;
                this.factors.SVD_V = V;
            elseif strcmpi(this.factorisation, 'chol')
                %figure;imagesc(K);axis image
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
            this.SF = SF;

            %this.lml = -GP_LMLG_FN(this, [this.meanpar, this.covpar, this.noisepar]);
            this.has_been_solved = 1;
        end


    % Query the model               
        function [mu_star, var_star, var_full] = query(this, x_star, qtask)

            if (~this.has_been_solved)
                error('GP must be solved first using GP.solve.');
            end
            if (~this.check())
                error('Abort.');
            end

            SF = this.SF;
            N = size(this.K,1); 
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
                warning('X is very large - querying data points in series.\n');
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

            ntask = this.ntask;
            mu_0 = this.MeanFn{qtask}.eval(x_star, this.meanpar);

            % build taskmap and nx
            nxp = 0;
            nt = this.ntask;
            for t = 1:nt
                nxp = nxp + size(this.X{t},2);
            end
            taskmap = zeros(1,nxp);
            upto = 0;
            for t = 1:nt
                nthis = size(this.X{t},2);
                taskmap(upto+(1:nthis)) = t;
                upto = upto + nthis;
            end

            for i = 1:step:nx % oops stopping too early - nx is really nq
                % progress bar
                if (large_query)
                    fprintf('%d...',i);
                end
                LR = i:min(nx, i+step-1);
                x_star_i = x_star(:,LR);
                ks = zeros(size(x_star_i,2),N);
                for t = 1:ntask
                    thistask = (taskmap ==t);
                    ks(:,thistask) = SF(qtask,t)*this.CovFn.eval(this.X{t},x_star_i,this.covpar{t}, this.covpar{qtask})';
                end

                product = (ks*this.alpha)';

                mu_star(1,LR) = mu_star(1,LR) + mu_0(1,LR) + product;

                if (nargout==2) % dont compute this if nargout = 1 (no variance) or 3 (full variance)
                    var0 = SF(qtask,qtask)*this.CovFn.pointval(x_star(:,LR), this.covpar{qtask});
                    if use_svd
                        % We got S2 from S2 = S2(:,ones(1,size(x_star(:,LR),2)));
                        v = bsxfun(@times, S2, (ks*this.factors.SVD_U)');
                    elseif use_chol
                        v = this.factors.L\ks';
                    end
                    predvar = max(0,var0 - sum(v.*v));
                    var_star(1,LR) = predvar;
                end
            end
            fprintf('\n');

            % if the user wants var_full then we need the whole covariance
            % matrix rather than jsut its diagonal.
            if (nargout ==3)
                if (large_query)
                    error('Full posterior covariance matrix is too large to compute!');
                end
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
        function par = learn(this)

            if (~this.check())
                error('Abort.');
            end
            par0 = cat(2,this.sigpar,this.meanpar{:}, this.covpar{:}, this.noisepar{:});
            
            % There is some bizzare matlab bug here - this.objectfun is a
            % member variable of type function pointer and not actually a method,
            % yet matlab makes no distinction and automatically passes
            % a pointer to 'this' as the first argument... might need
            % fixing in the future by changing this.objectfun(this,theta)
            par = this.solver_function(@(theta) this.objectfun(theta), par0, this.opts);

            % Parameters are automatically unpacked by calling an
            % optimisation criterion: calling LML also lets us update the
            % marginal likelihood indicator...
            this.lml = -GP_MT_LMLG_FN(this,par);
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
            f_star   = zeros(1,nx);
            [mu_star var_star var_full] = this.query(x_star);
            LFull = chol(var_full+eye(size(var_full))*1e-8, 'lower');
            fFull = [LFull*randn(1,nx)']';
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
            if (size(this.y,1) ~= 1)
                error('Y is transposed?\n');
            end
            valid = false;
            ntask = this.ntask;

            if ( size(this.sigpar,1)~=1)
                error('Sigpar is transposed?\n');
            end
            if (size(this.sigpar,2)~=ntask*(ntask+1)/2)
                error('Sigpar should be 1x(nt*(nt+1)/2)');
            end

            use_svd = strcmpi(this.factorisation, 'svd');
            use_chol = strcmpi(this.factorisation, 'chol');
            if ((use_svd==0)&&(use_chol==0))
                fprintf('Matrix factorisation should either be SVD or CHOL\n');
                drawnow; 
                return;
            end

            if ~isa(this.CovFn,'GP_MTCovFunc')
                error('Not a valid Covariance Function\n');
            end
            ncovpar = this.CovFn.npar(D); % Should be the same for each task

            for i=1:ntask
                if ~isa(this.MeanFn{i},'GP_MeanFunc')
                    error('Not a valid Mean Function Class\n');
                end
                if ~isa(this.NoiseFn{i},'GP_NoiseFunc')
                    error('Not a valid Noise Function Class\n');
                end
                if (size(this.covpar{i},1) ~= 1)
                    error('Each covpar should be 1*npar.\n');
                end 
                if (size(this.covpar{i},2) ~= ncovpar)
                    error('covpar is the wrong length.\n');
                end
                if (size(this.meanpar{i},1) > 1)
                    error('meanpar should be 1*npar.\n');
                end 
                nmeanpar = this.MeanFn{i}.npar(D);
                if (size(this.meanpar{i},2) ~= nmeanpar)
                    error('meanpar is the wrong length.\n');
                end    
            end

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
                params = cat(2,this.sigpar,this.meanpar{:}, this.covpar{:}, this.noisepar{:});
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
        end 
        
    end %methods
    
end % class   