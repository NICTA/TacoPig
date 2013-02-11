=======
TacoPig   
=======

TacoPig is a simple, object-oriented Gaussian process MATLAB toolbox.

Model Structure
===============

tacopig/TacoPigStructure.png illustrates the structure of the Gaussian process
model used by TacoPig. 


.. image:: https://github.com/NICTA/TacoPig/raw/master/TacoPigStructure.png

At the model's centre is a Gaussian process class instance which contains the
key functions for inference and learning as well as the properties such as the
training data. Typical examples of this type of class are the
tacopig.gp.Regressor and tacopig.gp.ClassifierLaplace (not ye implemented)
which perform regression and classification, respectively.

The GP class instance needs 4 plug-ins to operate.

1. Mean Function
2. Covariance Function
3. Noise Function
4. Objective Function

A Gaussian process is defined completely by a mean function and a covariance
function. A Gaussian Process can be written as :math:`f(x)~GP(m(x),k(x,x^*))`
Where the function :math:`f` is distributed as a Gaussian process with mean
function :math:`m(x)` and covariance function :math:`k(x,x*)`.  Noise is often
present in real data and needs to be accounted for to enable the model to learn
the true underlying function.  Learning the hyperparameters of the Gaussian
process is often done by minimising some objective function.

Getting Started 
===============

Add the TacoPig tool box to MATLAB's path.  Demos can be found in tacopig/demos
Open the demos and run them to see how the Gaussian process model is
initialised, trained and queried.


Creating a Gaussian Process Instance
===============
TacoPig consists of a number of classes that can be combined to create a Gaussian 
process instance in a similar manner to that illustrated in Figure 1. Figure 2 
shows the range of classes implemented in Version 0.1 of TacoPig.
 
.. image:: https://github.com/NICTA/TacoPig/raw/master/TacoPig_v0Point1.png

A Gaussian process instance can be initialised by combining at least one block from each of the 5 class groups (Gaussian Processes, Mean Functions, Covariance Functions, Noise Functions and Objective Functions)

For example, Figure 3 shows the structure of a Gaussian process regressor instantiation with a stationary mean and noise functions. Its covariance function is the squared exponential and during parameter/hyperparameter learning, the optimiser attempts to minimise the negative log marginal likelihood.
 
.. image:: https://github.com/NICTA/TacoPig/raw/master/TacoPig_Simple.png

The corresponding MATLAB code::

   GP = tacopig.gp.Regressor;
   % Plug in the components
   GP.MeanFn  = tacopig.meanfn.StationaryMean();
   GP.CovFn   = tacopig.covfn.SqExp();
   GP.NoiseFn = tacopig.noisefn.Stationary();
   GP.objective_function = @tacopig.objectivefn.NLML;


Figure 4 shows a more complicated structure of a Gaussian process regressor instantiation with a fixed mean function (its parameters will not change during the learning phase) and a LogStationary noise function with certain parameters clamped (in the code below, the first parameter of the LogStationary noise function is clamped to 1.2). The covariance function consists of a summing of the exponential covariance function and a Matern5 covariance function with some of its hyperparameters remapped (in the code below the first and second hyperparameters are remapped to always have the same value, the third hyperparameter is independent.). Cross-validation is used by the optimiser during the learning phase.

.. image:: https://github.com/NICTA/TacoPig/raw/master/TacoPig_Complex.png


The corresponding MATLAB code::

   GP = tacopig.gp.Regressor;
   % Plug in the components
   GP.MeanFn  = tacopig.meanfn.FixedMean();
   GP.CovFn   = tacopig.covfn.Sum(tacopig.covfn.Exp(),tacopig.covfn.Remap(tacopig.covfn.Mat5(),[1 1 2]));
   GP.NoiseFn = tacopig.noisefn.Clamp(tacopig.noisefn.LogStationary(),[1],[1.2]);
   GP.objective_function = @tacopig.objectivefn.CrossVal;


Code Documentation
==================

To view the code documentation, type:

``doc tacopig``

into the MATLAB command line after you have added the toolbox to MATLAB's path.


Contributing to TacoPig
=======================

Bugs
----

Please submit bug reports using the GitHub Issues feature. Be sure to first check
that the bug exists in the latest release, as it may have been fixed since you last
downloaded TacoPig. We strongly encourage you to include code to reproduce the
error if possible.

Contributing Code
-----------------

We welcome contributions to TacoPig in the form of GitHub pull requests. Before
you create a pull request, please ensure you have included the following:

1. Descriptive commit messages
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Describe what functionality you have added. If you have fixed a bug,
refer to the issue number using a hash: e.g. "fixed #453". The standard
format for a commit message is a short one-line description of the commit,
then a blank line, then a more detailed paragraph if necessary.

2. Documentation
~~~~~~~~~~~~~~~~
Comments should be integrate with the MATLAB doc and help commands. Use the 
comments included in the MATLAB scripts of TacoPig v0.1 as templates.

3. Test cases
~~~~~~~~~~~~~

The majority of TacoPig contributions will add covariance, mean or noise
function objects to the code. A standard testing suite
(``tacopig.tests.standardtestsXML``) is provided for testing these plugin
objects, and will produce an XML report if a filename is provided as the first
argument.  standardtestsXML will automatically detect new plugin objects
provided they are placed in the appropriate sub-package folders.

This standard test suite checks for common problems such as outputs that are
the wrong size, covariance functions that are not positive definite, or
gradients that do not agree with numerical differentiation.

However, if you are contributing a plugin with non-standard functionality, or a
new type of GP core object, then it will be neccessary to provide explicit
additional testing. Additional scripts should be placed in +tacopig/+tests.
 


4. Your name as a contributor in the AUTHORS file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Contact
=======

* Alistair: alistair.reid@nicta.com.au
* Simon: simon.ocallaghan@nicta.com.au
* Lachlan: lachlan.mccalman@nicta.com.au


