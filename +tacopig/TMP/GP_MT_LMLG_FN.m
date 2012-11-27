% Log marginal likelihood gradient for Multi-Task GPs
% (C) Alistair Reid & NICTA Sept 2012

function [nlml, nlmlg] = GP_MT_LMLG_FN(GP, parvec)
            
            % Special case: Multi-task
            if ~isa(GP,'GP_MultiTask')
                 error('Please use GP_LMLG_FN instead.');
            end
            % Get current solver configuration
            use_svd = strcmpi(GP.factorisation, 'svd');
            use_chol = strcmpi(GP.factorisation, 'chol');
            
            % Unpack the parameters
            % We actually update the GP hypers during learning in case we have to stop the optimisation early
            % Inverse of packing with cat(2,GP.sigpar,GP.meanpar{:}, GP.covpar{:}, GP.noisepar{:})
            ntask = GP.ntask;
            nsigs = ntask*(ntask+1)/2;
            GP.sigpar = parvec(1:nsigs);
            upto = nsigs;
            D = zeros(1,ntask);
            for i=1:ntask
                D(i) = size(GP.X{i},1);
            end
            for i=1:ntask
                nmpar = GP.MeanFn{i}.npar(D(i));
                GP.meanpar{i}  = parvec(upto + (1:nmpar));
                upto = upto + nmpar;
            end
            for i=1:ntask
                ncovpar = GP.CovFn.npar(D(i));
                GP.covpar{i}   = parvec(upto + (1:ncovpar));
                upto = upto + ncovpar;
            end
            for i=1:ntask
                nnoisepar = GP.NoiseFn{i}.npar(D(i));
                GP.noisepar{i} = parvec(upto + (1:nnoisepar));
                upto = upto + nnoisepar;
            end
            npars = upto; % count how many parameters we used in total
            
            % Apply these parameters to the model:
            L = zeros(ntask);
            L(tril(ones(ntask))>0.5) = GP.sigpar(:);
            SF = L*L';          
            X = GP.X;
            GP.SF = SF;
            
            % build taskmap and nx
            nx = 0;
            for t = 1:ntask
                nx = nx + size(X{t},2);
            end
            taskmap = zeros(1,nx);
            upto = 0;
            for t = 1:ntask
                nthis = size(X{t},2);
                taskmap(upto+(1:nthis)) = t;
                upto = upto + nthis;
            end
            K = zeros(nx); % Allocate Space
            Kbase = cell(ntask,ntask); % store the parts of K just in case
            mu = zeros(1,nx);
            covpars = GP.covpar;
            for i=1:ntask
                thistask = (taskmap==i);
                xi = X{i};
                mu(thistask) = GP.MeanFn{i}.eval_y(xi, GP.meanpar{i});
                pari = covpars{i};
                K0 = GP.CovFn.eval(xi, xi, pari, pari);
                Kbase{i,i} = K0;
                noise = GP.NoiseFn{i}.eval(xi, GP.noisepar{i});
                K(thistask,thistask) = SF(i,i)*K0 + noise;
                for j = i+1:ntask  
                    thattask = (taskmap==j);
                    xj = X{j};
                    parj = covpars{j};
                    Kbase{i,j} = GP.CovFn.eval(xi, xj, pari, parj); % storing for gradient
                    K0 = SF(i,j)*Kbase{i,j};
                    K(thistask,thattask) = K0;
                    K(thattask,thistask) = K0';
                end
            end   
            ym = (cat(2,GP.y{:}) - mu);
            GP.K = K;
            
            % Factorise K & compute -LML
            if (use_svd)
                [U,S,V] = svd(K);
                S0 = diag(S);
                S2 = S0;
                S2(S2>0) = 1./S2(S2>0);
                invK = V*diag(S2)*U';
                alpha = (invK*ym');
                nlml = 0.5*(ym*alpha + sum(log(S0)) + nx*log(2*pi));
            elseif (use_chol)
                L = chol(K, 'lower');
                if nargout>1
                    I = eye(size(K));
                    invK  = L'\(L\I);
                    alpha = (invK*ym');
                else
                    alpha = L'\(L\ym');
                end
                nlml = 0.5*(ym*alpha + 2*sum(log(diag(L))) + nx*log(2*pi));
            else
                error('Unrecognised factorisation choice');
            end
            
            
            % Compute the gradient for covariance hyperparameters
            % Much trickier than vanilla case due to the block-wise setup
            if nargout>1

                % Allocate a gradient vector...
                nlmlg = zeros(npars,1); % has to be tall to be compatible wit the framework

                % Cache the sensitivity of Kinv to K...
                W = invK - alpha*alpha'; 

                
                % In order of appearance in the parameter vector:
                
                % 1. NLML Gradient wrt sigvar hyperparameters
                L = tril(ones(ntask));
                [jj,ii] = meshgrid(1:ntask, 1:ntask);
                sampi = ii(L>0);
                sampj = jj(L>0);
                L(L>0) = GP.sigpar(:);
                % W=LL' - deriv(i,j) = jth column of L added to the ith 
                % row and its transpose added to the ith column of zeros(nt)
                for s = 1:nsigs
                    i = sampi(s);
                    j = sampj(s);
                    col = L(:,j);
                    deriv = zeros(ntask);
                    deriv(:,i) = col;
                    deriv(i,:) = col;
                    deriv(i,i) = deriv(i,i)*2;
                    gradmat = zeros(nx); % Gradient of K wrt hyper #S
                    
                    for i=1:ntask
                        thistask = (taskmap==i);
                        gradmat(thistask,thistask) = deriv(i,i)*Kbase{i,i};
                        for j = i+1:ntask
                            thattask = (taskmap==j);
                            K0 = Kbase{i,j};
                            gradmat(thistask,thattask) = deriv(i,j)*K0;
                            gradmat(thattask,thistask) = deriv(j,i)*K0';
                        end
                    end
                    % Compute NLMLG from kernel gradient
                    nlmlg(s) = 0.5*sum(sum(W.*gradmat)); 
                end
                upto = nsigs; % Count hypers computed dynamically

                
                % 2. NLML Gradient wrt meanfn...
                for i=1:ntask
                    thistask = (taskmap==i);
                    mgrad = GP.MeanFn{i}.gradient(X{i}, GP.meanpar{i});
                    nmpar = length(mgrad); %GP.MeanFn{i}.npar(D(i));
                    thisalpha = alpha(thistask);
                    for j=1:nmpar
                        nlmlg(upto + j) = -mgrad{j}*thisalpha; % dont compute over zeros
                    end
                    upto = upto + nmpar;
                end
                
                
                % 3. NLML wrt covfn hypers
                % This is easier to compute out of order
                
                % insert(task, hyper) gives the position in the gradient vector
                insert = cell(1,ntask);
                for i=1:ntask
                    nmpar = GP.CovFn.npar(D(i));
                    % use GP.meanpar{i}, insert{i} to reference
                    insert{i} = upto+(1:nmpar);
                    upto = upto + nmpar;
                end
                covpars = GP.covpar;
                for i=1:ntask
                    thistask = (taskmap==i);
                    xi = X{i};
                    pari = covpars{i};
                    for j = i:ntask  
                        thattask = (taskmap==j);
                        xj = X{j};
                        parj = covpars{j};
                        [g1,g2] = GP.CovFn.gradient(xi, xj, pari, parj);
                        
                        thisW1 = SF(i,j)*W(thistask,thattask);
                        thisW2 = SF(j,i)*W(thattask,thistask);
                        
                        % propagate these out
                        for k = 1:length(g1)
                            indx = insert{i}(k);   % which #hyper are we actually using?
                            nlmlg(indx) = nlmlg(indx) + 0.5*(sum(sum(thisW1.*g1{k})) + sum(sum(thisW2.*(g1{k}'))) );
                        end
                        
                        % Here is some madness:
                        % on the case i==j we end up double dipping on the
                        % two blocks thisW1 and thisW2 actually being the
                        % same block. However in this case it is also true
                        % that g1 = g2 due to them both being the same
                        % parameters.. therefore instead of having a
                        % special case for the {i,i} case, we can just not
                        % do the second gradient and it will halve the
                        % output exactly cancelling out our double dip on
                        % the above block :D obvious...
                        if (j>i)
                            for k = 1:length(g2)
                                indx = insert{j}(k);   % which #hyper are we actually using?
                                nlmlg(indx) = nlmlg(indx) + 0.5*(sum(sum(thisW1.*g2{k})) + sum(sum(thisW2.*(g2{k}'))) );
                            end
                        end
                    end
                end   
                
                % 4. NLML Gradient wrt noise parameters
                for i=1:ntask
                    ngrad = GP.NoiseFn{i}.gradient(X{i}, GP.noisepar{i});
                    nnpar = length(ngrad); %GP.MeanFn{i}.npar(D(i));
                    % noise gradient only applies to block diagonals
                    thistask = (taskmap==i);
                    thisW = W(thistask,thistask); 
                    for j=1:nnpar
                        nlmlg(upto+j) = 0.5*sum(sum(thisW.*ngrad{j}));
                    end
                    upto = upto + nnpar;
                end
                
                % 5. Do we spit out the gradients as we optimise?
                % This behaviour can slow things down and should default to off
                if (GP.spit_out_gradients)
                    GP.gcount = GP.gcount + 1;
                    if (GP.gcount>8)
                        GP.gcount = 0;
                        disp(parvec);
                        GP.check_gradients(parvec);
                    end
                end
                
            end
        end