%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Multi Task Gaussian Process Demo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% Generate some joint-Gaussian multi-task data
nx = 30;     % number of observations
nplot = 200; % number of points to plot
X0 = linspace(0,10, nplot);
xstar = X0;
r = randperm(rs, nplot);
r1 = sort(r(1:nx/2));
r2 = sort(r(nx/2+1:nx));
x1 = xstar(r1);
x2 = xstar(r2);

% Parameters of model:
ntask = 2;
LS1 = 0.50;         % Length Scale 1
LS2 = 0.45;         % Length Scale 2
SF = [1.0  0.9;     % Signal Variance
      0.9  1.0];
noise = 1e-3;       % Noise Standard Deviation
% Simulate some observations:
covfn = GP_MTSqExpCov();
K11  = covfn.eval(X0,X0,LS1,LS1);
K12  = covfn.eval(X0,X0,LS1,LS2);
K22  = covfn.eval(X0,X0,LS2,LS2);
K0   = [SF(1)*K11  SF(3)*K12; 
        SF(2)*K12' SF(4)*K22];
K0 = K0 + noise^2*eye(size(K0));
rand('state',1);
L = chol(K0,'lower');
samps = L*randn(rs,2*nplot,1);
Y1 = samps(1:nplot);
Y2 = samps(nplot+1:end);
y1 = Y1(r1)';
y2 = Y2(r2)';
figure
plot(X0,Y1,'r'); hold on
plot(x1,y1,'ro');
plot(X0,Y2,'b'); hold on
plot(x2,y2,'bo');
axis equal
title('Ground truth and Observations');


%% Define a multi-task GP model
GP          = GP_MultiTask(2);
GP.X        = {x1,x2};
GP.y        = {y1,y2};
GP.MeanFn   = {GP_ConstantMean(0), GP_ConstantMean(0)};  % Per Task Mean
GP.CovFn    = GP_MTSqExpCov();                           % Inter-Task Covariance
GP.NoiseFn  = {GP_StatNoise(), GP_StatNoise()};          % Per Task Noise
GP.optimset('gradobj', 'on');                           % No Gradients (yet)
GP.factorisation = 'SVD';

%% Seed initial hyperparameters with auto-sized vectors
ntask = 2;
ref  = chol(eye(ntask), 'lower');
GP.sigpar = ref(find(tril(ones(ntask))))';
GP.covpar   = cell(1,ntask);
GP.meanpar  = cell(1,ntask);
GP.noisepar = cell(1,ntask);
for i=1:ntask
    D = size(GP.X{i},1);
    GP.covpar{i}   = 0.25*ones(1, GP.CovFn.npar(D));
    GP.meanpar{i}  = 0.25*ones(1, GP.MeanFn{i}.npar(D));               % Constant Mean has no parameters
    GP.noisepar{i} = 0.25*ones(1, GP.NoiseFn{i}.npar(D)); 
end

% We could cheat to get these parameters...
% ref = chol(SF);
% GP.sigpar   = find(triu(ones(ntask)));
% GP.covpar   = {LS1,   LS2};
% GP.noisepar = {noise, noise};


%% Query the GP pre-training
GP.solve();

% to test gradients
% GP.optimset('gradobj','on');
% GP.check_gradients();
% return


[mf1, vf1] = GP.query(xstar,1);
[mf2, vf2] = GP.query(xstar,2);
mfn = {mf1,mf2};
vfn = {vf1,vf2};
cols = {'r','b'};
figure
hold on
sp = [0 0];
% note: we have 2 sigma bounds.. 
for i=1:2
    sp(i) = subplot(2,1,i);
    mf = mfn{i}; vf = vfn{i}; sf  = 2*sqrt(vfn{i});
    f  = [mf+2*sf,flipdim(mf-2*sf,2)]';
    h(1) = fill([xstar, flipdim(xstar,2)], f, [6 6 6]/8, 'EdgeColor', cols{i},'FaceAlpha', 0.5);
    hold on;
    if (i==1)
        h(2) = plot(X0,Y1,'k');
    else
        h(2) = plot(X0,Y2,'k');            
    end
    h(3) = plot(xstar,mf,[cols{i},'-'],'LineWidth',2);
    h(4) = plot(x1, y1, 'ro', 'MarkerSize', 5, 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor',[1 0 0]);
    h(5) = plot(x2, y2, 'bo', 'MarkerSize', 5, 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor',[0 0 1]);
    %legend(h,'Predictive Standard Deviation','True Function', 'Predictive Mean', 'Task 1 Training Points', 'Task 2 Training Points')
    axis equal
    title(sprintf('Prediction Task %d (not learnt)',i));
    axis tight
end
linkaxes(sp,'xy');

%% Learn & Re-Query
GP.learn();
GP.solve();
[mf1, vf1] = GP.query(xstar,1);
[mf2, vf2] = GP.query(xstar,2);
mfn = {mf1,mf2};
vfn = {vf1,vf2};
cols = {'r','b'};
figure
hold on
sp = [0 0];
for i=1:2
    sp(i) = subplot(2,1,i);
    mf = mfn{i}; vf = vfn{i}; sf  = 2*sqrt(vfn{i});
    f  = [mf+2*sf,flipdim(mf-2*sf,2)]';
    h(1) = fill([xstar, flipdim(xstar,2)], f, [6 6 6]/8, 'EdgeColor', cols{i},'FaceAlpha', 0.5);
    hold on;
    if (i==1)
        h(2) = plot(X0,Y1,'k');
    else
        h(2) = plot(X0,Y2,'k');            
    end
    h(3) = plot(xstar,mf,[cols{i},'-'],'LineWidth',2);
    h(4) = plot(x1, y1, 'ro', 'MarkerSize', 5, 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor',[1 0 0]);
    h(5) = plot(x2, y2, 'bo', 'MarkerSize', 5, 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor',[0 0 1]);
    %legend(h,'Predictive Standard Deviation','True Function', 'Predictive Mean', 'Task 1 Training Points', 'Task 2 Training Points')
    axis equal
    title(sprintf('Prediction Task %d (Learnt)',i));
    axis tight
end
linkaxes(sp,'xy');


% Give the gradients a good workout
if (false)
    for i=1:20
        GP.covpar = {0.5*rand,0.5*rand};GP.check_gradients();
    end
end