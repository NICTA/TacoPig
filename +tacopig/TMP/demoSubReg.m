%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Subset of Regressors Demo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initially, it is a seperate class, we can merge this functionality into
% other classes as needed... this is to explore the method

% Set up TPGP environment
clear all
clear functions
clc
close all
p = pwd;
root = strfind(p, 'geotherML')+9;
slash = p(root); root = p(1:root);
addpath(genpath([root, 'lib']));
rs = RandStream('mt19937ar','Seed',50); % make result repeatable


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

xstar = linspace(-8, 8, 201); % and query on a grid
[indxs, induced] = kmeans(X, 10); % only half as many points
% we will now compute the regression over these induced points

%% Configure the class model
GP = GP_SR();
GP.X = X;
GP.XI = induced'; % Not X, but clustered X to fewer points.
GP.y = y;
GP.factorisation = 'Chol'; 
GP.MeanFn =GP_ConstantMean(mean(y));
GP.CovFn = GP_SqExpCov();
%GP.NoiseFn = GP_StatNoise();
GP.optimset('gradobj', 'off');
% GP.objective_function should now default to the subset of regressors quasi-likelihood GP_SRLMLG_FN;


%% Seed initial hyperparameters with auto-sized vectors
GP.covpar  = 2*ones(1,GP.CovFn.npar(size(X,1)));
GP.meanpar = zeros(1,GP.MeanFn.npar(size(X,1)));
GP.noisepar = 0.4;% *ones(1,GP.NoiseFn.npar);


%% Solve and query
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
h(4) = plot(induced', interp1(xstar,mf,induced'), 'ro');
title('Not learnt');
legend(h,'Predictive Standard Deviation','Predictive Mean', 'Training Points', 'Induced Points')

%% Learn, solve & Query
disp('Press any key to begin learning.')
pause
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
title('Learnt');
legend(h,'Predictive Standard Deviation','Predictive Mean', 'Training Points', 'Induced Points')
