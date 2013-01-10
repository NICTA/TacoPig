%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Script to test regressor on the 1d example problem
%  Instead of visualising, it runs a test suite and generates a report
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
fprintf('Testing GP Regression framework___________________________________\n');
% Make sure we have the optimization toolbox available
addpath(genpath(['../optimization']))
report = JenkinsReport('Tacopig');
suite = JenkinsSuite('Test_Regressor');
report.addsuite(suite);

% 1-D Testing Data
% Rasmussen & Williams "Gaussian Processes for Machine Learning", Fig. 2.5
X = [-7.3005, -6.3772, -6.3481, -5.8868, -4.7591, -4.0911, -3.7461,...
    -2.7995, -2.1775, -0.9235, 0.7502, 1.0004, 2.3851, 2.4582,...
    4.2504, 4.3248, 4.8933, 5.8177, 6.1426];
y = [ -3.1593, -2.6634, -2.6364, -2.5123, -2.0282, -2.2231, -2.4161,...
    -2.1390, -1.5605, -1.4464, -0.4838, -0.2821, -1.4956,...
    -1.6061, -2.8353, -2.7971, -2.5217, -2.2265, -2.1399];
xstar = linspace(-8, 8, 51); 
[X id] = sort(X);
eqthresh = 1e-2; % equality threshold for approximately equal


%% Minimum Setup for the GP
% Use a standard GP regression model:
good_setup = true;
thistest = JenkinsTest('Regressor','Setup');
try
    GP = tacopig.gp.Regressor;
    GP.X = X;
    GP.y = y;
    GP.MeanFn  = tacopig.meanfn.ConstantMean(mean(y));
    GP.CovFn   = tacopig.covfn.SqExp();
    GP.NoiseFn = tacopig.noisefn.Stationary();
    GP.objective_function = @tacopig.objectivefn.NLML;
    GP.covpar   = 0.5*ones(1,GP.CovFn.npar(size(X,1)));
    GP.meanpar  = zeros(1,GP.MeanFn.npar(size(X,1)));
    GP.noisepar = 1e-1*ones(1,GP.NoiseFn.npar);
    GP.check();
    GP.opts.Display = 'off';
catch e
    str = geterr(e);
    thistest.fail = 1;
    thistest.message = sprintf('Failed initial setup.', current);
    thistest.description = str;
    good_setup = false;
end
suite.addtest(thistest)

% GP.objective_function = @tacopig.objectivefn.CrossVal;

 mftarget = ...
      [ -2.4421   -2.8265   -3.1016   -3.0602   -2.8307   -2.6582   -2.5660   -2.4449...
        -2.2723   -2.1157   -2.0334   -2.0567   -2.1915   -2.3616   -2.4431   -2.3887...
        -2.1979   -1.8898   -1.6161   -1.5332   -1.5522   -1.5093   -1.4661   -1.5624...
        -1.7164   -1.6750   -1.2716   -0.6640   -0.3401   -0.5922   -1.0969   -1.3983...
        -1.4965   -1.6320   -1.8497   -2.0647   -2.2936   -2.5722   -2.7807   -2.7661...
        -2.5773   -2.3912   -2.2905   -2.2304   -2.1582   -2.0898   -2.0579   -2.0552...
        -2.0592   -2.0612   -2.0617];

    vftarget = ...
       [0.2148    0.1117    0.0127    0.0454    0.0509    0.0057    0.0102    0.0190...
        0.0893    0.0745    0.0109    0.0211    0.0111    0.0070    0.0295    0.0543...
        0.0131    0.0208    0.0106    0.0537    0.1436    0.0973    0.0109    0.0753...
        0.1890    0.2001    0.1040    0.0182    0.0070    0.0521    0.1433    0.1282...
        0.0287    0.0192    0.1245    0.2122    0.2089    0.1136    0.0144    0.0137...
        0.0109    0.0301    0.0467    0.0118    0.0073    0.0466    0.1605    0.2334...
        0.2487    0.2500    0.2500];

% Pre-Learning Query
if (good_setup)
    GP.solve(); 
    [mf, vf] = GP.query(xstar);
   
    thistest = JenkinsTest('Regressor','Predictive mean pre-learning.');
    if (approxequal(mf,mftarget)>eqthresh)
        thistest.fail = 1;
        thistest.message = 'Regression Test Fail';
        thistest.description = 'GP output mean (pre learning) has changed from historical value.';
        fprintf('Problem: %s\n', thistest.description);
    end
    suite.addtest(thistest);
    
    thistest = JenkinsTest('Regressor','Predictive variance pre-learning.');
    if (approxequal(vf,vftarget)>eqthresh)
        thistest.fail = 1;
        thistest.message = 'Regression Test Fail';
        thistest.description = 'GP output variance (pre learning) has changed from historical value.';
        fprintf('Problem: %s\n', thistest.description);
    end
      
    suite.addtest(thistest);
    
else
    thistest = JenkinsTest('Regressor','Predictive mean pre-learning.');
    thistest.skip = 1;
    suite.addtest(thistest);
    thistest = JenkinsTest('Regressor','Predictive variance pre-learning.');
    thistest.skip = 1;
    suite.addtest(thistest);
end


% Pre-Learning Query
if (good_setup)
    GP.factorisation = 'SVD';
    GP.solve(); 
    [mf, vf] = GP.query(xstar);
   
    thistest = JenkinsTest('Regressor','Predictive mean pre-learning with SVD factorisation.');
    if (approxequal(mf,mftarget)>eqthresh)
        thistest.fail = 1;
        thistest.message = 'Regression Test Fail';
        thistest.description = 'GP output mean (pre learning) has changed from historical value.';
        fprintf('Problem: %s\n', thistest.description);
    end
    suite.addtest(thistest);
    
    thistest = JenkinsTest('Regressor','Predictive variance pre-learning with SVD factorisation.');
    if (approxequal(vf,vftarget)>eqthresh)
        thistest.fail = 1;
        thistest.message = 'Regression Test Fail';
        thistest.description = 'GP output variance (pre learning) has changed from historical value.';
        fprintf('Problem: %s\n', thistest.description);
    end
      
    suite.addtest(thistest);
    
else
    thistest = JenkinsTest('Regressor','Predictive mean pre-learning.');
    thistest.skip = 1;
    suite.addtest(thistest);
    thistest = JenkinsTest('Regressor','Predictive variance pre-learning.');
    thistest.skip = 1;
    suite.addtest(thistest);
end



thistest = JenkinsTest('Regressor','Learn with NLML and minfunc.');
if (good_setup)
    GP.factorisation = 'CHOL';
    GP.learn();
    GP.solve();
    err = approxequal([GP.covpar, GP.meanpar,GP.noisepar],[1.1119    0.7811  -0.0097]);
    if (err>eqthresh)
        thistest.fail = 1;
        thistest.message = 'Learning Reached different Optimum';
        thistest.description = 'GP learnt hyperparameters have changed from historical value.';
        fprintf('Problem: %s\n', thistest.description);       
    end
else
    thistest.skip = 1;
end
suite.addtest(thistest);
% Ensure the results of this test dont affect further tests
GP.covpar = [1.1119  0.7811];
GP.meanpar = [];
GP.noisepar = -0.0097;
GP.solve();

%% Generate samples from prior and posterior

thistest = JenkinsTest('Regressor','Sample prior with fixed seed.');
if (good_setup)
    fstar = GP.sampleprior(xstar, 50);
    target = ...
        [-2.0297   -2.0026   -1.8119   -1.5518   -1.2218   -0.9516   -0.8255   -0.9456...
         -1.3268   -1.8303   -2.2912   -2.5351   -2.4657   -2.2386   -2.0068   -1.8487...
         -1.7708   -1.7526   -1.7698   -1.8249   -1.9917   -2.2354   -2.4936   -2.7630...
         -2.9708   -3.1180   -3.1808   -3.1143   -2.9801   -2.8042   -2.6496   -2.4687...
         -2.2653   -1.9842   -1.6764   -1.3456   -1.0826   -0.8258   -0.6612   -0.6313...
         -0.8008   -1.2337   -1.8788   -2.6062   -3.2318   -3.6609   -3.8443   -3.7907...
         -3.5030   -3.0512   -2.4514];

     if (approxequal(fstar,target)>eqthresh)
        thistest.fail = 1;
        thistest.message = 'Prior sample has unexpected realisation.';
        thistest.description = 'Prior sample has changed from historical value.';
        fprintf('Problem: %s\n', thistest.description);       
    end 
else
    thistest.skip = 1;
end
suite.addtest(thistest);

thistest = JenkinsTest('Regressor','Sample posterior with fixed seed.');
if (good_setup)
    fstar = GP.sampleposterior(xstar, 50);
    target = ...
        [-3.3448   -3.3287   -3.1678   -2.9875   -2.7723   -2.6610   -2.5796   -2.4957   -2.3547   -2.1563   -2.0280   -2.0285   -2.1607   -2.3902   -2.5087   -2.4279   -2.1966   -1.9036...
         -1.6212   -1.4226   -1.3724   -1.3864   -1.4377   -1.5112   -1.4786   -1.3090   -0.9856   -0.5933   -0.3206   -0.2412   -0.4201   -0.7638   -1.2526   -1.7597   -2.2419   -2.5897...
         -2.8252   -2.8755   -2.8625   -2.7466   -2.5878   -2.4558   -2.3459   -2.2536   -2.1564   -2.0910   -2.0797   -2.1123   -2.1007   -2.0137   -1.7543];
    if (approxequal(fstar,target)>eqthresh)
        thistest.fail = 1;
        thistest.message = 'Posterior sample has unexpected realisation.';
        thistest.description = 'Posterior sample has changed from historical value.';
        fprintf('Problem: %s\n', thistest.description);       
    end 
else
    thistest.skip = 1;
end
suite.addtest(thistest);

thistest = JenkinsTest('Regressor','Check LML Function.');
if (good_setup)
    theta = [GP.meanpar, GP.covpar, GP.noisepar];
    [LML,LMLg] = GP.objectfun(theta);

    if (approxequal(LML,-2.6014)>eqthresh)
        thistest.fail = 1;
        thistest.message = 'Objective function problem.';
        thistest.description = 'Objective function has changed from historical value.';
        fprintf('Problem: %s\n', thistest.description);       
    elseif (numel(LML)~=1)
        thistest.fail = 1;
        thistest.message = 'Objective function problem.';
        thistest.description = 'Objective function should be scalar when called without gradients.';
        fprintf('Problem: %s\n', thistest.description);       
    end 
else
    thistest.skip = 1;
end
suite.addtest(thistest);

thistest = JenkinsTest('Regressor','Check LML gradient.');
if (good_setup)
    if (approxequal(LMLg,[-0.0271; 0.0050; -2.3410])>eqthresh)
        thistest.fail = 1;
        thistest.message = 'Objective function gradient problem.';
        thistest.description = 'Objective function has changed from historical value.';
        fprintf('Problem: %s\n', thistest.description); 
    end
else
   thistest.skip = 1;
end 
fprintf('__________________________________________________________________\n');
   


    % Checks we know it does for sure but might want to test anyway:
    % bad factorisation
    % MeanFn is a tacopig.meanfn.MeanFunc
    % CovFn is a tacopig.covfn.CovFunc
    % covpar is 1xn and the correct size
    % correct parameter lengths
