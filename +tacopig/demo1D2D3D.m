%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Gaussian Process Demo Script
%  Demonstrates GP regression using the taco-pig toolbox.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1-D Example%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
close all; clear all; clear functions; clc;
% import tacopig.*;

%% Set up 1-D Data
% Training Data
% Rasmussen & Williams "Gaussian Processes for Machine Learning", Fig. 2.5
X = [-2.1775,-0.9235,0.7502,-5.8868,-2.7995,4.2504,2.4582,6.1426,...
    -4.0911,-6.3481,1.0004,-4.7591,0.4715,4.8933,4.3248,-3.7461,...
    -7.3005,5.8177,2.3851,-6.3772];
y = [1.4121,1.6936,-0.7444,0.2493,0.3978,-1.2755,-2.221,-0.8452,...
    -1.2232,0.0105,-1.0258,-0.8207,-0.1462,-1.5637,-1.098,-1.1721,...
    -1.7554,-1.0712,-2.6937,-0.0329];
xstar = linspace(-8, 8, 201); 

% Order data points for visualisation
[X id] = sort(X);
y = y(id);

%% Set up Gaussian process

% Use a standard GP regression model:
GP = tacopig.gp.Regressor;

% Plug in the data
GP.X = X;
GP.y = y;

% Plug in the components
GP.MeanFn  = tacopig.meanfn.ConstantMean(mean(y));
% GP.CovFn   = tacopig.covfn.Sum(tacopig.covfn.Mat3(),tacopig.covfn.SqExp());%SqExp();
GP.CovFn   = tacopig.covfn.SqExp();
GP.NoiseFn = tacopig.noisefn.Stationary();
GP.objective_function = @tacopig.objectivefn.NLML;
GP.solver_function = @anneal;

% Initialise the hyperparameters
GP.covpar   = 0.5*ones(1,GP.CovFn.npar(size(X,1)));
GP.meanpar  = zeros(1,GP.MeanFn.npar(size(X,1)));
GP.noisepar = 1e-3*ones(1,GP.NoiseFn.npar);


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
title('Not learnt');
legend(h,'Predictive Standard Deviation','Predictive Mean', 'Training Points')

%% Learn & Query
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
title('Learnt');
legend(h,'Predictive Standard Deviation','Predictive Mean', 'Training Points')

%% Generate samples from prior and posterior
figure; subplot(1,2,1)
xstar = linspace(-8,8,100);
hold on;
for i = 1:5
    fstar = GP.sampleprior(xstar);
    plot(xstar,fstar, 'color', rand(1,3));
end
title('Samples from Prior')

subplot(1,2,2)
plot(X, y, 'k+', 'MarkerSize', 17)
xstar = linspace(-8,8,100);
hold on;
for i = 1:5
    fstar = GP.sampleposterior(xstar);
    plot(xstar,fstar, 'color', rand(1,3));
end
title('Samples from Posterior')
pause


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2-D Example%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
close all; clear all; clear functions; clc;
% import tacopig.*;

%% Set up 2-D Data
% Training Data
X = (rand(2,30)-0.5)*6;
y = peaks(X(1,:),X(2,:))+1e-2*randn(1,size(X,2));

[xeva yeva] = meshgrid(-3:0.2:3,-3:0.2:3);
xstar = [xeva(:)';yeva(:)'];

figure; scatter3(X(1,:),X(2,:),y,40,y,'filled')


%% Set up Gaussian process

% Use a standard GP regression model:
GP = tacopig.gp.Regressor;

% Plug in the data
GP.X = X;
GP.y = y;

% Plug in the components
GP.MeanFn  = tacopig.meanfn.ConstantMean(mean(y));
GP.CovFn   = tacopig.covfn.SqExp();%SqExp();
GP.NoiseFn = tacopig.noisefn.Stationary();
GP.objective_function = @tacopig.objectivefn.NLML;
GP.solver_function = @anneal;

% Initialise the hyperparameters
GP.covpar   = 1*ones(1,GP.CovFn.npar(size(X,1)));
GP.meanpar  = zeros(1,GP.MeanFn.npar(size(X,1)));
GP.noisepar = 1e-3*ones(1,GP.NoiseFn.npar);


%% Learn & Query
GP.learn();
GP.solve();
[mf, vf] = GP.query(xstar);
sf  = sqrt(vf);

% Display learnt model
figure
scatter3(X(1,:),X(2,:),y,60,y,'filled')
hold on
surf(xeva,yeva,reshape(mf,size(xeva)))
title('Predictive Mean Function and Standard Deviation Surfaces');
hold on
surf(xeva,yeva,reshape(mf+sf,size(xeva)),'facealpha',0.1)
surf(xeva,yeva,reshape(mf-sf,size(xeva)),'facealpha',0.1)
pause

%% Generate samples from prior and posterior
figure; 
hold on;
for i = 1:5
    clf
    fstar = GP.sampleprior(xstar);
    surf(xeva,yeva,reshape(fstar,size(xeva)));
    title('Samples from Prior')
    pause(1)
end


figure
hold on;
for i = 1:5
    clf
    fstar = GP.sampleposterior(xstar);
    surf(xeva,yeva,reshape(fstar,size(xeva)));
    title('Samples from Posterior')
    pause(1)
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% 3-D Example%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
close all; clear all; clear functions; clc;
% import tacopig.*;

%% Set up 3-D Data
% Training Data
groundtruth = @(x,y,z) 5*exp(-(((x).^2)/5+((1-y).^2)/2+((0.5-z).^2)/3))...
    -4*exp(-(((2-x).^2)/2+((-1-y).^2)/5+((-1.5-z).^2)/2));
X = (rand(3,100)-0.5)*6;
y = groundtruth(X(1,:),X(2,:),X(3,:))+1e-2*randn(1,size(X,2));

[xeva yeva zeva] = meshgrid(-3:0.2:3,-3:0.2:3,-3:0.2:3);
xstar = [xeva(:)';yeva(:)';zeva(:)'];

figure; scatter3(X(1,:),X(2,:),X(3,:),40,y,'filled')


%% Set up Gaussian process

% Use a standard GP regression model:
GP = tacopig.gp.Regressor;

% Plug in the data
GP.X = X;
GP.y = y;

% Plug in the components
GP.MeanFn  = tacopig.meanfn.ConstantMean(0);
GP.CovFn   = tacopig.covfn.SqExp();%SqExp();
GP.NoiseFn = tacopig.noisefn.Stationary();
GP.objective_function = @tacopig.objectivefn.NLML;
GP.solver_function = @anneal;

% Initialise the hyperparameters
GP.covpar   = 1*ones(1,GP.CovFn.npar(size(X,1)));
GP.meanpar  = zeros(1,GP.MeanFn.npar(size(X,1)));
GP.noisepar = 1e-3*ones(1,GP.NoiseFn.npar);


%% Learn & Query
GP.learn();
GP.solve();
[mf, vf] = GP.query(xstar);
sf  = sqrt(vf);

%% Visualise GP Outputs

% Build a colormap
cmap = jet(5); 
levels = linspace(-2.5,3.5,5);

% Generate isosurfaces
figure; 
for ii = 1:5
    camlight
    lighting gouraud
    hh(ii) = patch(isosurface(xeva, yeva, zeva, reshape(mf,size(xeva)), levels(ii)));
    set(hh(ii), 'Facecolor', cmap(ii,:), 'Edgecolor', 'none', 'facealpha', (1-(5-abs(levels(ii)))/5));
    axis([-3 3 -3 3 -3 3])
    caxis([-4,5]);
    xlabel('x');ylabel('y');zlabel('z')
    colorbar
end
 title('Predictive Mean')






% % Display learnt model
% figure
% scatter3(X(1,:),X(2,:),y,60,y,'filled')
% hold on
% surf(xeva,yeva,reshape(mf,size(xeva)))
% title('Predictive Mean Function and Standard Deviation Surfaces');
% hold on
% surf(xeva,yeva,reshape(mf+sf,size(xeva)),'facealpha',0.1)
% surf(xeva,yeva,reshape(mf-sf,size(xeva)),'facealpha',0.1)
% pause
% 
% %% Generate samples from prior and posterior
% figure; 
% hold on;
% for i = 1:5
%     clf
%     fstar = GP.sampleprior(xstar);
%     surf(xeva,yeva,reshape(fstar,size(xeva)));
%     title('Samples from Prior')
%     pause(1)
% end
% 
% 
% figure
% hold on;
% for i = 1:5
%     clf
%     fstar = GP.sampleposterior(xstar);
%     surf(xeva,yeva,reshape(fstar,size(xeva)));
%     title('Samples from Posterior')
%     pause(1)
% end

