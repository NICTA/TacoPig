%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                       Gaussian Process Demo
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear all; clear functions; clc;
% Add working path
p = pwd;
root = strfind(p, 'geotherML')+9;
slash = p(root);     root = p(1:root);
Libroot = [root, 'lib', slash];
addpath(genpath(Libroot))

% Use the same data as Rasmussen in Figure 2.5
X = [-2.1775;-0.9235;0.7502;-5.8868;-2.7995;4.2504;2.4582;6.1426;...
    -4.0911;-6.3481;1.0004;-4.7591;0.4715;4.8933;4.3248;-3.7461;...
    -7.3005;5.8177;2.3851;-6.3772]';
y = [1.4121;1.6936;-0.7444;0.2493;0.3978;-1.2755;-2.221;-0.8452;...
    -1.2232;0.0105;-1.0258;-0.8207;-0.1462;-1.5637;-1.098;-1.1721;...
    -1.7554;-1.0712;-2.6937;-0.0329]';
n = size(X,2);
[X id] = sort(X);
y = y(id);
% y = y+X; % Adds a linear trend if uncommented
% y(X>=0) = y(X>=0)+5;
xstar = linspace(-8, 8, 201); % and query on a grid


%% Configure the standard GP model
GP = GP_STD(); 
GP.X = X;
GP.y = y;
GP.factorisation = 'SVD'; % 'SVD' or 'CHOL' or 'JITCHOL'

%% Mean Functions:
% Some examples of mean functions
%       GP_LinearMean();
%       GP_ConstantMean( mean(y)); 
    GP.MeanFn =GP_ConstantMean(mean(y));

%% Covariance Functions:
% List of Covariance functions
%       GP_SqExpCov(); 
%       GP_ExpCov(); 
%       GP_Mat3Cov(); 
%       GP_Mat5Cov();
%       GP_CovSum(covfn1,covfn2,covfn3...); 
%       GP_CovProd(covfn1,covfn2,covfn3...); 
%       GP_Projection(S, covfn); 
%       GP_ClampCov(covfn, params_index, params_values); 
%       GP_Remap(covfn, remap_vector);  

% Some examples of non-trivial covariance functions
% Sum of two covariance functions: GP_CovSum(GP_SqExpCov(), GP_ExpCov());
% Product of two covariance functions: GP_CovProd(GP_SqExpCov(), GP_ExpCov());
% Project the observations through a sensitivity matrix (inverse problem): GP_Projection(Grav_Sensitivity, GP_SqExpCov())
% Clamp the second hyperparameter to a value of 1.3: GP_ClampCov(GP_ExpCov(),2,1.3)
% Both hypers index the same optimisation parameter: GP_Remap(GP_SqExpCov(), [1 1])
% GP.CovFn = GP_Remap(GP_CovSum(GP_ExpCov(),GP_SqExpCov()), [1 1 2 2]); % learn a sqexp with length scale = variance

%     GP.CovFn = GP_SqExpCov_Barrier();
        GP.CovFn = GP_SqExpCov();

    
%% Noise Functions:
% List of noise functions
%       GP_StatNoise(); % Standard stationary iid noise
%       GP_LinIncreaseNoise(); %Noise linearly increases with abs(X(1,:))
%       GP_ClampNoise(noisefn, params_index, params_values); 

    GP.NoiseFn = GP_StatNoise();


%% Optimisation Parameters
% Choose optimiser % defaults to fminunc
    GP.solver_function = @fminunc; 
% Choose Objective function  % defaults to gplmlgrad
     GP.objective_function = @GP_LMLG_FN;
  %  GP.objective_function = @GP_CRSVAL_FN; % cross validation
% Use analytical gradients?
    GP.optimset('gradobj', 'off');


%% Seed initial hyperparameters with auto-sized vectors
GP.covpar  = 0.5*ones(1,GP.CovFn.npar(size(X,1)));
GP.meanpar = zeros(1,GP.MeanFn.npar(size(X,1)));
GP.noisepar = 1e-3*ones(1,GP.NoiseFn.npar);
%% Solve before learning
GP.solve(); 
[mf, vf] = GP.query(xstar);
sf  = sqrt(vf);

% Display before learning
figure
plot(X, y, 'k+', 'MarkerSize', 17)
f  = [mf+2*sf,flipdim(mf-2*sf,2)]';
h(1) = fill([xstar, flipdim(xstar,2)], f, [6 6 6]/8, 'EdgeColor', [6 6 6]/8);
hold on
h(2) = plot(xstar,mf,'k-','LineWidth',2);
h(3) = plot(X, y, 'k+', 'MarkerSize', 17);
title('Not learnt');
legend(h,'Predictive Standard Deviation','Predictive Mean', 'Training Points')

%% Learn & Query
disp('Press any key to begin learning.')
pause
GP.learn();
GP.solve();
% GP.K(end/2+1:end,1:end/2-1)=0
% GP.K(1:end/2-1,end/2+1:end)=0
[mf, vf] = GP.query(xstar);
sf  = sqrt(vf);


% Display learnt model
figure
plot(X, y, 'k+', 'MarkerSize', 17)
f  = [mf+2*(sf),flipdim(mf-2*(sf),2)]';
h(1) = fill([xstar, flipdim(xstar,2)], f, [6 6 6]/8, 'EdgeColor', [6 6 6]/8);
hold on
h(2) = plot(xstar,mf,'k-','LineWidth',2);
h(3) = plot(X, y, 'k+', 'MarkerSize', 17);
title('Learnt');
legend(h,'Predictive Standard Deviation','Predictive Mean', 'Training Points')

%% Generate samples from prior and posterior
figure; subplot(1,2,1)
xstar = linspace(-8,8,100);
hold on;
for i = 1:10
    fstar = GP.sampleprior(xstar);
    plot(xstar,fstar, 'color', rand(1,3));
end
title('Samples from Prior')

subplot(1,2,2)
plot(X, y, 'k+', 'MarkerSize', 17)
xstar = linspace(-8,8,100);
hold on;
for i = 1:10
    fstar = GP.sampleposterior(xstar);
    plot(xstar,fstar, 'color', rand(1,3));
end
title('Samples from Posterior')

% Output GP object to console
GP

% figure;
% for i = 1:100
%     BarrierLoc = (i-50)*8/50
%    GP.covpar(end-1) = BarrierLoc 
%     GP.solve();
%     lml(i) = GP.lml
% end
% figure;plot(lml)