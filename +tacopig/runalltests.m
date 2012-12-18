function runalltests(varargin)
    % make sure we change to the right directory...
    mfilepath = mfilename('fullpath');
    mfilepath = mfilepath(1:end-length(mfilename())-length('+TacoPig/'));
    testdir = [mfilepath,'matlab-xunit-doctest/xunit']
    addpath(testdir)
    runtests('tacopig.tests')
    return;
