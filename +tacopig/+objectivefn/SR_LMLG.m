% Gaussian Process LMLG function for subset of regressors
% note: there is is a matlab bug causing this to behave as a class member
% function...
function [nlml, nlmlg] = SR_LMLG(this, parvec)
    % a dummy gradient
    nlml = zeros(size(parvec));

    % unpack the hyperparameter vector
    D = size(this.X,1);
    ncovpar = this.CovFn.npar(D);
    nmeanpar = this.MeanFn.npar(D);
    nnoisepar = 1;
    
    this.meanpar = parvec(1:nmeanpar);
    this.covpar = parvec(nmeanpar+1:nmeanpar+ncovpar);
    this.noisepar = parvec(nmeanpar+ncovpar+1:nmeanpar+ncovpar+nnoisepar);
    
    % cut and pasted from solve
    mu = this.MeanFn.eval_y(this.X, this.meanpar);
    KI  = this.CovFn.Keval(this.XI, this.covpar);
    KIX = this.CovFn.eval(this.XI, this.X, this.covpar);

    % add a tiny noise to KI to keep positive definiteness
    eps = 1e-6*sum(diag(KI)); % or could use min etc
    KI  = KI + eps*eye(size(KI));

    N = size(this.X,2);  % number of training points
    m = size(this.XI,2); % number of induced points
    noise = this.noisepar;
    
    ym = (this.y - mu);
    
    
%     We could do it this way, except that figuring out pseudoK and then
%     inverting is just as slow as using all the regressors... so we do
%     some numerical tricks below....
%     pseudoK = KIX'*KI*KIX + noise^2*eye(N) ;
%     L = chol(pseudoK); 
%     nlml = 0.5*ym*(L'\(L\ym')) + sum(log(diag(L))) + 0.5*N*log(2*pi);
%     if(nargout==2)
%         error('Gradients not implemented');
%         nlmlg = cell(1,parvec);
%         for i=1:parvec
%             nlmlg{i} = 0;
%         end
%     end
    
    
    L     = chol(KI,'lower');
    V     = L\KIX;
    
    
    pK    = V*V'+(noise^2)*eye(m); %pseudoK
    Lm    = chol(pK,'lower');
    b     = (Lm\V)*ym';

    % Compute the log marginal likelihood
    yb    = ym*ym'-b'*b;

    nlml  = sum(log(diag(Lm))) + 0.5*(N-m)*log(noise^2) + ...
            0.5*(1/(noise^2))*yb + 0.5*N*log(2*pi);
    

    
    if(nargout==2)
        error('Gradients not implemented');
        nlmlg = cell(1,parvec);
        for i=1:parvec
            nlmlg{i} = 0;
        end
    end
    
    return
    
%     
%     % if we need to get gradient
%         alpha = -2*yt' + 2*(c'*c)*yt';
%     %Compute necessary matrices
%         Lt    = L*Lm;
%         B1    = Lt'\(Lm\V);
%         %Compute the variance
%         %opts.TRANSA = true;
%         %opts.LT     = true;
%         %v           = linsolve(L,ks,opts);
% 
%         %B1 = linsolve(Lt,linsolve(Lm,V,opts),opts);
%         %B1 = solve_tril(Lt,solve_tril(Lm,V));
%         b1    = Lt'\b;
%         invKI = L'\(L\I);
%         invA  = Lt'\(Lt\I);
% 
%         %Compute the gradient
%         KId   = feval(kgfun,kpar,KI,X(:,si));
%         tmp   = feval(kgfun,kpar,KIn,X(:,si),X);
%         mu    = V'*(Lm'\b);
%         ng    = zeros(1,length(KId)+1);
% 
%         for i = 1:length(tmp)
%             ng(i) = -0.5*sum(sum((invKI-(noise^2)*invA).*KId{i})) + ...
%                     sum(sum(B1.*tmp{i})) - (noise^(-2))*b1'*tmp{i}*(yt-mu')'...
%                     + 0.5*b1'*KId{i}*b1;
%         end
%         %Noise gradient
%         ng(end) = ((noise^2)*trace(Lm'\(Lm\I)) + n - m + ...
%                   sum((Lm'\b).^2) - (noise^(-2))*yb)*(noise)^(-1);