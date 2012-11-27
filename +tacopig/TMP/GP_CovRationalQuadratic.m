par = log(par);

ell = exp(par(1));
sf2 = exp(2*par(2));
alpha = exp(par(3));

% precompute squared distances
if dg                                                               % vector kxx
  D2 = zeros(size(x,1),1);
else
  if xeqz                                                 % symmetric matrix Kxx
    D2 = sq_dist(x'/ell);
  else                                                   % cross covariances Kxz
    D2 = sq_dist(x'/ell,z'/ell);
  end
end

if nargin<4                                                        % covariances
  K = sf2*((1+0.5*D2/alpha).^(-alpha));
else                                                               % derivatives
  if i==1                                               % length scale parameter
    K = sf2*(1+0.5*D2/alpha).^(-alpha-1).*D2;
  elseif i==2                                              % magnitude parameter
    K = 2*sf2*((1+0.5*D2/alpha).^(-alpha));
  elseif i==3
    K = (1+0.5*D2/alpha);
    K = sf2*K.^(-alpha).*(0.5*D2./K - alpha*log(K));
  else
    error('Unknown hyperparameter')
  end
end