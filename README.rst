=======
TacoPig   
=======

TacoPig is a simple, object-oriented Gaussian process MATLAB toolbox.

Model Structure
===============

tacopig/TacoPigStructure.png illustrates the structure of the Gaussian process
model used by TacoPig. 


.. image:: https://github.com/NICTA/TacoPig/raw/develop/TacoPigStructure.png

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
*insert paragraph on documentation here*

3. Test cases
~~~~~~~~~~~~~
*insert paragraph on test cases here*

4. Your name as a contributor in the AUTHORS file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Contact
=======

* Alistair: alistair.reid@nicta.com.au
* Simon: simon.ocallaghan@nicta.com.au

