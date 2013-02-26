%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Gaussian Process Demo Script
%  Demonstrates GP regression using the taco-pig toolbox on 1-D Data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Add optimization folder
if ~exist('minfunc')
    fprintf('It looks like you need to add minfunc to your default path...\n');
    fprintf('(Add tacopig/optimization/minfunc{/mex} to pathdef.m for permanent access)\n');
    fprintf('Press any key to attempt to continue...\n');
    pause();
end
tacopigroot = which('tacopig.taco');
tacopigroot = tacopigroot(1:end-15);
addpath(genpath([tacopigroot,'optimization']))

%% %%%%%%%%%%%%%% 1-D Example of Subset of Regressors%%%%%%%%%%%%%%%%%%%%% 
close all; clear functions; clc;

%% Set up 1-D Data
% Training Data
% Rasmussen & Williams "Gaussian Processes for Machine Learning", Fig. 2.5
X = [-2.1775,-0.9235,0.7502,-5.8868,-2.7995,4.2504,2.4582,6.1426,...
    -4.0911,-6.3481,1.0004,-4.7591,0.4715,4.8933,4.3248,-3.7461,...
    -7.3005,5.8177,2.3851,-6.3772];
y = [1.4121,1.6936,-0.7444,0.2493,0.3978,-1.2755,-2.221,-0.8452,...
    -1.2232,0.0105,-1.0258,-0.8207,-0.1462,-1.5637,-1.098,-1.1721,...
    -1.7554,-1.0712,-2.6937,-0.0329];

n = size(X,2);
[X id] = sort(X);
y = y(id);

xstar = linspace(-8, 8, 201); 
try % pick the induced points. only half as many points in this case.
    [indxs, induced] = kmeans(X, n/2);
catch
%     induced = (rand(10,1)-0.5) *abs(max(X)-min(X))+mean(X) ;
    induced = [linspace(min(X),max(X),n/2)]';
end
% we will now compute the regression over these induced points


%% Set up Gaussian process

% Use a standard GP regression model:
GP = tacopig.gp.SubsetRegressor;

% Plug in the data
GP.X = X;
GP.XI = induced';
GP.y = y;

% Plug in the components
GP.MeanFn  = tacopig.meanfn.FixedMean(mean(y));
GP.CovFn   = tacopig.covfn.SqExp();
GP.NoiseFn = tacopig.noisefn.Stationary();
GP.objective_function = @tacopig.objectivefn.SR_LMLG;

% GP.solver_function = @anneal;

% Initialise the hyperparameters
GP.covpar   = 2*ones(1,GP.CovFn.npar(size(X,1)));
GP.meanpar  = zeros(1,GP.MeanFn.npar(size(X,1)));
GP.noisepar = 0.4*ones(1,GP.NoiseFn.npar);


%% Before Learning: Query
GP.solve(); 
[mf, vf] = GP.query(xstar);
sf  = sqrt(vf);

% Display predicitve mean and variance
figure
plot(X, y, 'k+', 'MarkerSize', 17)
f  = [mf+2*sf,flipdim(mf-2*sf,2)]';
h(1) = fill([xstar, flipdim(xstar,2)], f, [6 6 6]/8, 'EdgeColor', [6 6 6]/8);
hold on
h(2) = plot(xstar,mf,'k-','LineWidth',2);
h(3) = plot(X, y, 'k+', 'MarkerSize', 17);
h(4) = plot(induced', interp1(xstar,mf,induced'), 'ro');
title('Before Hyperparameter Training');
legend(h,'Predictive Standard Deviation','Predictive Mean', 'Training Points','Induced Points','Location','SouthWest')

%% Learn & Query
GP.learn();
GP.solve();
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
h(4) = plot(induced', interp1(xstar,mf,induced'), 'ro');
title('After Hyperparameter Training');
legend(h,'Predictive Standard Deviation','Predictive Mean', 'Training Points','Induced Points','Location','SouthWest')
