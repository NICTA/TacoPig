% TacoPig Toolbox
% Version 0.1 (R2012a) 14-Dec-2012
%
%  Covariance Functions:
%   tacopig.covfn.Clamp             - Function that fixes specified hyperparameters of a covariance function during training
%   tacopig.covfn.CovFunc           - Covariance function abstract class
%   tacopig.covfn.Mat3              - Matern3 covariance function
%   tacopig.covfn.Mat5              - Matern5 covariance function
%   tacopig.covfn.Exp               - Exponential covariance function
%   tacopig.covfn.Product           - Covariance function that multiply two other covariance functions
%   tacopig.covfn.Remap             - Matern5 covariance function
%   tacopig.covfn.SqExp             - Squared exponential covariance function
%   tacopig.covfn.Sum               - Covariance function that sums two other covariance functions
%
%  Gaussian Core Mechanics:
%   tacopig.gp.GpCore               - Core properties of the Gaussian process
%   tacopig.gp.Regressor            - Gaussian process regression equations
%   
%  Mean Functions:
%   tacopig.meanfn.ConstantMean     - A constant value is used for the GP's mean
%   tacopig.meanfn.MeanFunc         - Mean function abstract class
%   tacopig.meanfn.StationaryMean   - A constant global mean value that is learnt during the hyperparameter training phase
%
%  Noise Functions
%   tacopig.noisefn.Clamp           - Function that fixes specified hyperparameters of a noise function during training
%   tacopig.noisefn.NoiseFunc       - Noise function abstract class
%   tacopig.noisefn.Stationary      - A constant global noise value that is learnt during the hyperparameter training phase
%   tacopig.noisefn.LogStationary   - Same as Stationary except the parameter that is passed to the noise function is the log of the true noise.
%
%  Objective Functions: Functions that are minimised during the training phase
%   tacopig.objectivefn.CrossVal    - k-fold cross validation (Currently set to k = 5)
%   tacopig.objectivefn.NLML        - Negative log marginal likelihood
%
% Last Revision: n/a.
% Copyright BSD Licence



