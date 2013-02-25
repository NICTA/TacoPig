% Gaussian Process LMLG function for subset of regressors
function [nlml, nlmlg] = SR_LMLG(this, parvec)
    % unpack the hyperparameter vector
    D = size(this.X,1);
    ncovpar = this.CovFn.npar(D);
    nmeanpar = this.MeanFn.npar(D);
    nnoisepar = this.NoiseFn.npar;
    
    this.meanpar = parvec(1:nmeanpar);
    this.covpar = parvec(nmeanpar+1:nmeanpar+ncovpar);
    this.noisepar = parvec(nmeanpar+ncovpar+1:nmeanpar+ncovpar+nnoisepar);
    
    % cut and pasted from solve
    mu = this.MeanFn.eval(this.X, this);
    KI  = this.CovFn.Keval(this.XI, this);
    KIX = this.CovFn.eval(this.XI, this.X, this);

    % add a tiny noise to KI to keep positive definiteness
    eps = 1e-6*sum(diag(KI)); % or could use min etc
    KI  = KI + eps*eye(size(KI));

    N = size(this.X,2);  % number of training points
    m = size(this.XI,2); % number of induced points
    noise = this.noisepar;
    
    ym = (this.y - mu);
    
    
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