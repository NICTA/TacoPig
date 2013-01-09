
% A standard test on all existing covariance functions
% If you have any specific tests, then you should write functions for them
% using Matlab xUnit.

% Standard testing with three suites:
% Detects all covariance functions and tests the following:
%   *  parse the code with matlab
%   *  instanciate the class with either no arguments or using the code
%      in the constant string member .ExampleUsage
%   *  Check that the npar method returns a value - uses this to make a
%      random test problem with a fixed seed reset so its repeatable
%
%   *  Check for a 'tacopig:dimMismatch' when X1 and X2 have different
%      dimensionality in eval(X1,X2,par)
%   *  Check the size of the output matrix of eval
%   *  Check that eval enforces the correct number of hyperparameters
%   *  Check that Keval enforces the correct number of hyperparameters
%   *  Check that pointval enforces the correct number of hyperparameters
%   *  Check that Keval returns a matrix of the correct size
%   *  Regression test: ensure Keval returns the same matrix as last time
%        - if there is no reference, then save automatically
%        - if it doesn't match, then ask the user whether we overwrite
%   *  Check size of eval(X,X,pars)
%   *  Check eval(X,X,pars) is same as Keval(X,pars)
%   *  Check that Keval(X,pars) is positive semidefinite (using chol)
%   *  Check pointeval(X,pars) is the diag of Keval(X,pars)
%   *  Warn if gradients are not functioning
%   *  If gradients are functioning:
%         - numerically check gradients at 10 random operating points
% Detects all mean functions and tests the following:
%   *  parse the code with matlab
%   *  instanciate the class with either no arguments or using the code
%      in the constant string member .ExampleUsage
%   *  Check that the npar method returns a value - uses this to make a
%      random test problem with a fixed seed reset so its repeatable
%   *  Check for a 'tacopig:dimMismatch' when X1 and X2 have different
%      dimensionality in eval(X,par)
%   *  Check that eval enforces the correct number of hyperparameters
%   *  Check the size of the output matrix of eval
%   *  Regression test: ensure eval returns the same matrix as last time
%        - if there is no reference, then save automatically
%        - if it doesn't match, then ask the user whether we overwrite
%   *  Warn if gradients are not functioning
%   *  If gradients are functioning:
%         - numerically check gradients at 10 random operating points
% Detects all noise functions and tests the following:
%   *  parse the code with matlab
%   *  instanciate the class with either no arguments or using the code
%      in the constant string member .ExampleUsage
%   *  Check that the npar method returns a value - uses this to make a
%      random test problem with a fixed seed reset so its repeatable
%   *  Check the size of the output matrix of eval
%   *  Regression test: ensure eval returns the same matrix as last time
%        - if there is no reference, then save automatically
%        - if it doesn't match, then ask the user whether we overwrite
%   *  Warn if gradients are not functioning
%   *  If gradients are functioning:
%         - numerically check gradients at 10 random operating points
% And builds a Jenkins-compatible XML report of the results.


clear all
clear functions
clc
testpath = '../+tacopig'; % Relative path of main library to here
interactive = false; % Do we prompt on regression test failures?
report = JenkinsReport('Tacopig');



%% Covariance function testing

fprintf('Testing Covariance Functions______________________________________\n\n');


covSuite = JenkinsSuite('Covariance Functions');
report.addsuite(covSuite);

% Find all the covariance functions
wd = pwd;
if wd(1) == '/'
    slash = '/';
else
    slash = '\';
end
cov_functions = [testpath, slash, '+covfn',slash,'*.m'];
filez = dir(cov_functions);
ncovs = length(filez)-1;
cov_names = cell(1, ncovs);
offset = 0;
for i=1:ncovs
    fname = filez(i+offset).name(1:end-2);
    if strcmp(fname,'CovFunc')
        offset = 1;
        fname = filez(i+offset).name(1:end-2);
    end
    cov_names{i} = fname;
end

% Compare to previous runs (regression test)
try
    load regression_test_cache.mat; % regressions
catch e
    fprintf('No regression tests found >> Creating new cache...\n');
    regressions = [];
end
eqthresh = 1e-4; % threshold for being approximately equal


for i=1:ncovs
    
    % Test 1 Create instances of these classes
    current = cov_names{i};
    thistest = JenkinsTest(current,'InitClass');
    
    
    fprintf('Testing %s:', current);
    pass = true;
    good_init = true;
    argcount = 0;
    try
        eval(['argcount = nargin(@tacopig.covfn.', current, ');']);
    catch e
        str = geterr(e);
        thistest.fail = 1;
        thistest.message = sprintf('Cannot initialise %s.', current);
        thistest.description = str;
        fprintf('\n%s\n', str);
        good_init = false;
    end
    
    currentFull = ['tacopig.covfn.', current];
    
    if (good_init)
        try
            if (argcount~=0)
                eval(['ExampleUsage = ',currentFull,'.ExampleUsage;']);
                if (isempty(ExampleUsage))
                    thistest.description = 'No test string provided';
                    error('failblog');
                end
                fullExampleUsage = ['covobj = ',ExampleUsage,';'];%tacopig.covfn.
                eval(fullExampleUsage);
            else
                eval(['covobj = ', currentFull,'();']);
            end
            
            rand('seed', 50);
            npoint = 15;
            ndims = 3;
            X  = rand(ndims,npoint);
            X2 = rand(ndims,npoint-1);
            X3 = rand(ndims-1,npoint);
            npar = covobj.npar(ndims);
            pars = rand(1,npar);
            pars2 = rand(1,npar-1);
            
        catch e
            fprintf('\nAborting Tests:%s', e.message);
            fprintf('\n%s\n',fullExampleUsage);
            %fprintf('Cannot continue testing %s:',current);
            thistest.fail = 1;
            thistest.message = sprintf('Cannot initialise %s.', current);
            pass = false;
            %if isempty(thistest.description)
            %thistest.description = geterr(e);
            %fprintf('\n%s >> %s\n',current,e.message);
            %end
            good_init = false;
        end
    end
    covSuite.addtest(thistest)
    
    % Look for standard exceptions
    thistest = JenkinsTest(current,'CheckEvalXInputDims');
    if (~good_init)
        thistest.skip = 1;
    else
        if ~throws(covobj, 'k = in.eval(varargin{1}, varargin{2}, varargin{3});', 'tacopig:dimMismatch', X, X3, pars)
            if (pass)
                fprintf('\n');
            end
            fprintf('%s.eval is not detecting input dimensionality mismatch between x1 and x2 and should throw tacopig:dimMismatch\n', current);
            pass = false;
            thistest.fail = 1;
            thistest.message = 'Assertion Failed.';
            thistest.description = [current,'.eval is not detecting input dimensionality mismatch between x1 and x2 and should throw tacopig:dimMismatch'];
        end
    end
    covSuite.addtest(thistest);
    
    % checking number of hyperparameters
    thistest = JenkinsTest(current,'CheckEvalVerifiesHyperLength');
    if (~good_init)
        thistest.skip = 1;
    else
        if ~throws(covobj, 'k = in.eval(varargin{1}, varargin{2}, varargin{3});', 'tacopig:inputInvalidLength', X, X2, [pars,1])
            if (pass)
                fprintf('\n');
            end
            fprintf('%s.eval is not detecting the wrong number of hypers and should throw tacopig:inputInvalidLength\n', current);
            thistest.fail = 1;
            thistest.message = 'Assertion Failed.';
            thistest.description = [current,'.Eval is not detecting the wrong number of hypers and should throw tacopig:inputInvalidLength.'];
            pass = false;
        end
    end
    covSuite.addtest(thistest);
    
    % checking number of hyperparameters on Keval
    thistest = JenkinsTest(current,'CheckKEvalVerifiesHyperLength');
    if (~good_init)
        thistest.skip = 1;
    else
        if ~throws(covobj, 'k = in.Keval(varargin{1}, varargin{2});', 'tacopig:inputInvalidLength', X, [pars,1])
            if (pass)
                fprintf('\n');
            end
            fprintf('%s.Keval is not detecting the wrong number of hypers and should throw tacopig:inputInvalidLength\n', current);
            thistest.fail = 1;
            thistest.message = 'Assertion Failed.';
            thistest.description = [current,'.Keval is not detecting the wrong number of hypers and should throw tacopig:inputInvalidLength.'];
            pass = false;
        end
    end
    covSuite.addtest(thistest);
    
    % checking number of hyperparameters on pointval
    thistest = JenkinsTest(current,'CheckpointvalVerifiesHyperLength');
    if (~good_init)
        thistest.skip = 1;
    else
        if ~throws(covobj, 'k = in.pointval(varargin{1}, varargin{2});', 'tacopig:inputInvalidLength', X, [pars,1])
            if (pass)
                fprintf('\n');
            end
            fprintf('%s.pointval is not detecting the wrong number of hypers and should throw tacopig:inputInvalidLength\n', current);
            thistest.fail = 1;
            thistest.message = 'Assertion Failed.';
            thistest.description = [current,'.pointval is not detecting the wrong number of hypers and should throw tacopig:inputInvalidLength.'];
            pass = false;
        end
    end
    covSuite.addtest(thistest);
    
    % check that Keval is the correct size
    thistest = JenkinsTest(current,'CheckKevalOutputDims');
    if (~good_init)
        thistest.skip = 1;
    else
        try
            K1 = covobj.Keval(X,pars);
            s = size(K1);
            if ~((s(1)==npoint)&(s(2)==npoint))
                error(sprintf('%s.Keval returned a matrix of the wrong size!',current));
            end
        catch e
            if (pass)
                fprintf('\n');
            end
            thistest.fail = 1;
            thistest.message = 'Output Wrong Size';
            thistest.description = sprintf('%s >> %s\n', current, e.message);
            fprintf('%s >> %s\n', current, e.message);
            pass = false;
        end
        covSuite.addtest(thistest);
        
        % and gives the same result as before
        thistest = JenkinsTest(current,'RegressionKevalOutput');
        if (~good_init)
            thistest.skip = 1;
        else
            try
                eval(sprintf('oldresult = regressions.cov.%s.Keval;',current));
            catch e
                fprintf('No previous result for %s.Keval - saving for next time...\n', current);
                eval(sprintf('regressions.cov.%s.Keval = K1;',current));
            end
            
            try
                if (approxequal(K1,oldresult)>eqthresh)
                    if (pass)
                        fprintf('\n');
                    end
                    fprintf('K1 doesnt match the cached regression data!\n');
                    
                    
                    if (interactive)
                        reply = input('Are you sure it is correct NOW? Y/N [Y]: ', 's');
                        if isempty(reply)
                            reply = 'N';
                        end
                        if strcmpi(reply,'y')
                            eval(sprintf('regressions.cov.%s.eval = M1;',current));
                            pass = true;
                        else
                            pass = false;
                        end
                        
                    end
                    
                    if (pass == false)
                        thistest.fail = 1;
                        thistest.message = 'Regression Test Fail';
                        thistest.description = 'Keval output does not match cached records.';
                    end
                end
            catch e
                thistest.fail = 1;
                thistest.message = 'Regression Test Fail';
                thistest.description = e.message;
            end
        end
    end
    covSuite.addtest(thistest);
    
    thistest = JenkinsTest(current,'KevalOutputPosDef');
    if (~good_init)
        thistest.skip = 1;
    else
        try
            % add a tiny representative noise
            K1 = K1 + 1e-4*diag(diag(K1));
            L = chol(K1);
        catch e
            if (pass)
                fprintf('\n');
            end
            fprintf('%s >> may not be positive definite!\n', current);
            pass = false;
            thistest.fail = 1;
            thistest.message = 'Regression Test Fail';
            thistest.description = 'Keval output is not positive definite!';
        end
    end
    covSuite.addtest(thistest);
    
    thistest = JenkinsTest(current,'KevalcompareEval');
    if (~good_init)
        thistest.skip = 1;
    else
        try
            K2 = covobj.eval(X,X,pars);
            Kpp = K2;
            s = size(K2);
            if ~((s(1)==npoint)&(s(2)==npoint))
                error('%s.eval returned a matrix of the wrong size!',current);
            end
            
            % compare K1 and K2 through covariance?
            if (approxequal(K1,K2)>eqthresh)
                error(sprintf('%s.Keval(X,par) and %s.eval(X,X,par) are not giving the same answer!\n',current, current));
            end
        catch e
            if (pass)
                fprintf('\n');
                pass = false;
            end
            fprintf('%s\n',e.message);
            thistest.fail = 1;
            thistest.message = 'Assertion failed.';
            thistest.description = 'Keval is not giving the same answer as eval.';
        end
    end
    covSuite.addtest(thistest);
    
    % check pointeval = diag
    thistest = JenkinsTest(current,'ComparePointevalWithEval');
    if (~good_init)
        thistest.skip = 1;
    else
        try
            target = diag(Kpp)';
            obtained = covobj.pointval(X,pars);
            if (approxequal(target,obtained) > eqthresh)
                error(sprintf('%s.pointeval(X,par) and %s.eval(X,X,par) are not giving consistent answers!\n',current, current));
            end
        catch e
            if (pass)
                fprintf('\n');
                pass = false;
            end
            thistest.fail = 1;
            thistest.message = 'Assertion failed.';
            thistest.description = e.message;
            fprintf('%s\n',e.message);
        end
        
    end
    covSuite.addtest(thistest);
    
    % Chack gradients if available - otherwise skip
    if (good_init)
        testgrads = true;
        try
            g = covobj.gradient(X, pars);
        catch e
            if (pass)
                fprintf('\n');
                pass = false;
            end
            fprintf('Gradients not functioning: %s\n',e.message);
            testgrads = false;
        end
    else
        testgrads = false;
    end
    
    thistest = JenkinsTest(current,'CheckGradientsAreCorrect');
    if ~(good_init&&testgrads)|(npar==0)
        thistest.skip = 1;
    else
        
        % If has gradients then check that they actually correspond to
        % working gradients
        
        try
            % check size of output vector
            if (size(g,1)~=1)
                error('Gradients should be 1xnpar!');
            end
            if (size(g,2)~=npar)
                error('Gradients should be 1xnpar!');
            end
            
            % Check that analytic is similar to numerical
            eps = 1e-6; % numerical step
            ntestpars = 10;
            errlist = '';
            for p = 1:ntestpars
                params = rand(1,npar);
                tX = rand(size(X));
                
                g = covobj.gradient(tX, params);
                
                %val0 = covobj.Keval(X, params);
                err = false;
                for i = 1:npar
                    params2 = params;
                    params3 = params;
                    params2(i) = params2(i) + eps;
                    params3(i) = params3(i) - eps;
                    val1 = covobj.Keval(tX,params2);
                    val2 = covobj.Keval(tX,params3);
                    numerical = (val1-val2)/(2*eps);
                    analytic  = g{i};
                    if (approxequal(target,obtained) > eqthresh)
                        fprintf('Gradient mismatch on parameter %d\n',i);
                        errlist = strcat(errlist, ' ', num2str(i));
                        err = true;
                    end
                end
                
            end
            if (err)
                error(['Gradients not matching numerical checks:',errlist]);
            end
            
        catch e
            if (pass)
                fprintf('\n');
                pass = false;
            end
            fprintf('%s\n',e.message);
            
            thistest.fail = 1;
            thistest.message = 'Assertion failed.';
            thistest.description = e.message;
        end
    end
    covSuite.addtest(thistest);
    
    if (pass)
        fprintf('Pass\n');
    else
        fprintf('\n')
    end
end
clear covobj;
fprintf('__________________________________________________________________\n\n');


% Mean Function Testing:
fprintf('\nTesting Mean Functions____________________________________________\n\n');
meanSuite = JenkinsSuite('Mean Functions');
report.addsuite(meanSuite);


% Find all the Mean Functions
mean_functions = [testpath, slash, '+meanfn',slash,'*.m'];
filez = dir(mean_functions);
nmf = length(filez)-1;
mfn_names = cell(1, nmf);
offset = 0;
for i=1:nmf
    fname = filez(i+offset).name(1:end-2);
    if strcmp(fname,'MeanFunc')
        offset = 1;
        fname = filez(i+offset).name(1:end-2);
    end
    mfn_names{i} = fname;
end

for i=1:nmf
    
    % Create instances of these classes
    current = mfn_names{i};
    thistest = JenkinsTest(current,'InitClass');
    fprintf('Testing %s:', current);
    pass = true;
    good_init = true;    
    argcount = 0;
    try
        eval(['argcount = nargin(@tacopig.meanfn.', current, ');']);
    catch e
        str = geterr(e);
        thistest.fail = 1;
        thistest.message = sprintf('Cannot initialise %s.', current);
        thistest.description = str;
        good_init = false;
        fprintf('\n%s\n', str);
        
    end
    
    currentFull = ['tacopig.meanfn.', current];
    
    if (good_init)
         try
            if (argcount~=0)
                eval(['ExampleUsage = ',currentFull,'.ExampleUsage;']);
                if (isempty(ExampleUsage))
                    thistest.description = 'No test string provided';
                    error('failblog');
                end
                fullExampleUsage = ['meanobj = ',ExampleUsage,';'];%tacopig.meanfn.
                eval(fullExampleUsage);
            else
                eval(['meanobj = ', currentFull,'();']);
            end
            
            rand('seed', 50);
            npoint = 15;
            ndims = 3;
            X  = rand(ndims,npoint);
            X2 = rand(ndims,npoint-1);
            X3 = rand(ndims-1,npoint);
            npar = meanobj.npar(ndims);
            pars = rand(1,npar);
            pars2 = rand(1,npar-1);
            
        catch e
            fprintf('\nAborting Tests:%s', e.message);
            fprintf('\n%s\n',fullExampleUsage);
            %fprintf('Cannot continue testing %s:',current);
            thistest.fail = 1;
            thistest.message = sprintf('Cannot initialise %s.', current);
            pass = false;
            %if isempty(thistest.description)
            %thistest.description = geterr(e);
            %fprintf('\n%s >> %s\n',current,e.message);
            %end
            good_init = false;
        end
    end
    meanSuite.addtest(thistest)

    
    % checking number of hyperparameters
    thistest = JenkinsTest(current,'CheckEvalVerifiesHyperLength');
    if (~good_init)
        thistest.skip = 1;
    else
        if ~throws(meanobj, 'k = in.eval(varargin{1}, varargin{2});', 'tacopig:inputInvalidLength', X, [pars,1])
            if (pass)
                fprintf('\n');
            end
            fprintf('%s.eval is not detecting the wrong number of hypers and should throw tacopig:inputInvalidLength\n', current);
            thistest.fail = 1;
            thistest.message = 'Assertion Failed.';
            thistest.description = [current,'.Eval is not detecting the wrong number of hypers and should throw tacopig:inputInvalidLength.'];
            pass = false;
        end
    end
    meanSuite.addtest(thistest); 
    
    % check that Keval is the correct size
    thistest = JenkinsTest(current,'CheckKevalOutputDims');
    if (~good_init)
        thistest.skip = 1;
    else
    
        % check that outputs are of the correct size
        try
            M1 = meanobj.eval(X,pars);
            s = size(M1);
            if ~((s(1)==1)&&(s(2)==npoint))
                error('%s.eval is the wrong size (should be 1*N)!',current);
            end
        catch e
            if (pass)
                fprintf('\n');
            end
            thistest.message = 'Assertion Failed';
            thistest.description = [current, '.eval is giving the wrong output dimensions.'];
            thistest.fail = 1;
            fprintf('%s >> %s\n', current, e.message);
            pass = false;
        end
    end
    meanSuite.addtest(thistest); 
        
        
    % and give the same result as before
    thistest = JenkinsTest(current,'RegressionEvalOutput');
    if (~good_init)
        thistest.skip = 1;
    else
        try
            eval(sprintf('oldresult = regressions.means.%s.eval;',current));
            if (approxequal(M1,oldresult)>eqthresh)
                if (pass)
                    fprintf('\n');
                end
                fprintf('eval doesn''t match the cached regression data!\n');

                reply = input('Are you sure it is correct NOW? Y/N [Y]: ', 's');
                if isempty(reply)
                    reply = 'N';
                end
                if strcmpi(reply,'y')
                    eval(sprintf('regressions.means.%s.eval = M1;',current));
                else
                    pass = false;
                end
                
                if (pass == false)
                    thistest.fail = 1;
                    thistest.message = 'Regression Test Fail';
                    thistest.description = 'Keval output does not match cached records.';
                end
                
            end
        catch e
            fprintf('No previous result for %s.eeval - saving for next time...\n', current);
            eval(sprintf('regressions.means.%s.eval = M1;',current));
        end
    end
    meanSuite.addtest(thistest);
        
    % Check if gradients are functioning
     if (good_init)
        testgrads = true;
        try
            g = meanobj.gradient(X, pars);
        catch e
            if (pass)
                fprintf('\n');
                pass = false;
            end
            fprintf('Gradients not functioning: %s\n',e.message);
            testgrads = false;
        end
    else
        testgrads = false;
    end
    
    
    thistest = JenkinsTest(current,'CheckGradientsAreCorrect');
    if ~(good_init&&testgrads)|(npar==0)
        thistest.skip = 1;
    else
    
        try
            g = meanobj.gradient(X, pars);
        catch e
            if (pass)
                fprintf('\n');
                pass = false;
            end
            fprintf('%s''s Gradients are not functioning: %s\n',current, e.message);
            testgrads = false;
        end

        % If has gradients then check that they actually correspond to working gradients
        if (testgrads&(npar>0))

            try
                % check size of output vector

                if ((size(g,1)~=1)||(size(g,2)~=npar))
                    error('Gradients should be 1xnpar!');
                end

                % Check that analytic is similar to numerical
                eps = 1e-6; % numerical step
                ntestpars = 10;

                for p = 1:ntestpars
                    params = rand(1,npar);
                    tX = rand(size(X));

                    g = meanobj.gradient(tX, params);

                    %val0 = meanobj.Keval(X, params);
                    err = false;
                    for i = 1:npar
                        params2 = params;
                        params3 = params;
                        params2(i) = params2(i) + eps;
                        params3(i) = params3(i) - eps;
                        val1 = meanobj.eval(tX,params2);
                        val2 = meanobj.eval(tX,params3);
                        numerical = (val1-val2)/(2*eps);
                        analytic  = g{i};
                        if (approxequal(numerical,analytic) > eqthresh)
                            fprintf('Gradient mismatch on parameter %d\n',i);
                            err = true;
                        end
                    end
                    if (err)
                        error('Possible gradient problem!');
                    end
                end

            catch e
                if (pass)
                    fprintf('\n');
                    pass = false;
                end
                fprintf('%s\n',e.message);
                thistest.fail = 1;
                thistest.message = 'Assertion failed.';
                thistest.description = e.message;
            end
        end
    end
    meanSuite.addtest(thistest);
    
    % passed everything?
    if (pass)
        fprintf('Pass\n');
    else
        fprintf('\n');
    end
end
clear meanobj;
fprintf('__________________________________________________________________\n\n');


% Find all the noise functions
fprintf('\nTesting Noise Functions:__________________________________________\n\n');
noiseSuite = JenkinsSuite('Noise Functions');
report.addsuite(noiseSuite);
noise_functions = [testpath, slash, '+noisefn',slash,'*.m'];
filez = dir(noise_functions);
nnf = length(filez)-1;
nfn_names = cell(1, nnf);
offset = 0;
for i=1:nnf
    fname = filez(i+offset).name(1:end-2);
    if strcmp(fname,'NoiseFunc')
        offset = 1;
        fname = filez(i+offset).name(1:end-2);
    end
    nfn_names{i} = fname;
end


for i=1:nnf
    
    % Create instances of these classes
    current = nfn_names{i};
    fprintf('Testing %s:', current);
    pass = true;
    
    
    thistest = JenkinsTest(current,'InitClass');
    good_init = true;
    
    argcount = 0;
    try
        eval(['argcount = nargin(@tacopig.noisefn.', current, ');']);
    catch e
        str = geterr(e);
        fprintf('\n%s\n', str);
        fprintf('Cannot continue testing %s\n',current);
        thistest.fail = 1;
        thistest.message = sprintf('Cannot initialise %s.', current);
        thistest.description = str;
        good_init = false;
    end
    currentFull = ['tacopig.noisefn.', current];
    
    
    if (good_init)
        try
            if (argcount~=0)
                eval(['ExampleUsage = ',currentFull,'.ExampleUsage;']);
                if (isempty(ExampleUsage))
                    thistest.description = 'No test string provided';
                    error('failblog');
                end
                fullExampleUsage = ['noiseobj = ',ExampleUsage,';'];%tacopig.covfn.
                eval(fullExampleUsage);
            else
                eval(['noiseobj = ', currentFull,'();']);
            end
            
            rand('seed', 50);
            npoint = 15;
            ndims = 3;
            X  = rand(ndims,npoint);
            X2 = rand(ndims,npoint-1);
            X3 = rand(ndims-1,npoint);
            npar = noiseobj.npar(ndims);
            pars = rand(1,npar);
            pars2 = rand(1,npar-1);
            
        catch e
            fprintf('\nAborting Tests:%s', e.message);
            fprintf('\n%s\n',fullExampleUsage);
            thistest.fail = 1;
            thistest.message = sprintf('Cannot initialise %s.', current);
            pass = false;
            good_init = false;
        end
    end
    noiseSuite.addtest(thistest)
    
    
    % Output sizes (could detect an accidental scalar output etc)
    thistest = JenkinsTest(current,'CheckevalOutputDims');
    if (~good_init)
        thistest.skip = 1;
    else
        try
            N1 = noiseobj.eval(X,pars);
            s = size(N1);
            if ~((s(1)==npoint)&&(s(2)==npoint))
                error('%s.eval is the wrong size (should be N*N)!',current);
            end
        catch e
            if (pass)
                fprintf('\n');
            end
            fprintf('%s >> %s\n', current, e.message);
            pass = false;
        end
    end
    noiseSuite.addtest(thistest);  
    
    % And regression test the output: 
    thistest = JenkinsTest(current,'RegressionKevalOutput');
    if (~good_init)
        thistest.skip = 1;
    else
        try
            eval(sprintf('oldresult = regressions.noise.%s.eval;',current));

            if (approxequal(N1,oldresult)>eqthresh)
                if (pass)
                    fprintf('\n');
                end
                fprintf('eval doesnt match the cached regression data!\n');

                reply = input('Are you sure it is correct NOW? Y/N [Y]: ', 's');
                if isempty(reply)
                    reply = 'N';
                end
                if strcmpi(reply,'y')
                    eval(sprintf('regressions.noise.%s.eval = N1;',current));
                else
                    pass = false;
                end
                
                 if (pass == false)
                    thistest.fail = 1;
                    thistest.message = 'Regression Test Fail';
                    thistest.description = 'eval output does not match cached records.';
                end
                
            end
        catch e
            fprintf('No previous result for %s.eval - saving for next time', current);
            eval(sprintf('regressions.noise.%s.eval = N1;',current));
        end
    end
    noiseSuite.addtest(thistest);
    
    
     % Chack gradients if available - otherwise skip
    if (good_init)
        testgrads = true;
        try
            g = noiseobj.gradient(X, pars);
        catch e
            if (pass)
                fprintf('\n');
                pass = false;
            end
            fprintf('Gradients not functioning: %s\n',e.message);
            testgrads = false;
        end
    else
        testgrads = false;
    end
    
    % If no gradients then warn
    thistest = JenkinsTest(current,'CheckGradientsAreCorrect');
    if ~(good_init&&testgrads)|(npar==0)
        thistest.skip = 1;
    else
        try
            % check size of output vector
            if ((size(g,1)~=1)||(size(g,2)~=npar))
                error('Gradients should be 1xnpar!');
            end
            
            % Check that analytic is similar to numerical
            eps = 1e-6; % numerical step
            ntestpars = 4;
            
            for p = 1:ntestpars
                params = rand(1,npar);
                tX = rand(size(X));
                
                g = noiseobj.gradient(tX, params);
                
                %val0 = noiseobj.Keval(X, params);
                err = false;
                for i = 1:npar
                    params2 = params;
                    params3 = params;
                    params2(i) = params2(i) + eps;
                    params3(i) = params3(i) - eps;
                    val1 = noiseobj.eval(tX,params2);
                    val2 = noiseobj.eval(tX,params3);
                    numerical = (val1-val2)/(2*eps);
                    analytic  = g{i};
                    if (approxequal(numerical,analytic) > eqthresh)
                        fprintf('Gradient mismatch on parameter %d\n',i);
                        err = true;
                    end
                end
                if (err)
                    error('Possible gradient problem!');
                end
            end
            
        catch e
            if (pass)
                fprintf('\n');
                pass = false;
            end
            fprintf('%s\n',e.message);
            thistest.fail = 1;
            thistest.message = 'Assertion failed.';
            thistest.description = e.message;
        end
    end
    noiseSuite.addtest(thistest);
    
    % passed everything?
    if (pass)
        fprintf('Pass\n');
    end
end
fprintf('__________________________________________________________________\n');

% Save the regression tests results as they stand
try
    save regression_test_cache.mat regressions; % regressions
catch e
    fprintf('Could not save regression test values... try saving manually!');
end

