% This function provides a standard testing suite
% New contributions to the covariance functions, mean functions or noise functions
% will be tested automatically for the basics below. Other contributions 
% may need new code for testing.
% The usage is given by:
%   tacopig.tests.standardtestsXML('filename.xml')
% If no filename is provided then no report is produced.
% The standard tests are applied to each object type:
% Covariance functions:
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
% Mean funcitons:
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
% Noise Functions:
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

function standardtestsXML(reportpath)


    testpath = mfilename('fullpath');
    % now cut off the +tests
    ind = findstr(testpath, '+tests');
    ind = ind(end);
    testpath = testpath(1:ind-2);

    % Find all the covariance functions
    wd = pwd;
    if wd(1) == '/'
        slash = '/';
    else
        slash = '\';
    end

    % sets up the report and everything
    interactive = false; % Do we prompt on regression test failures?
    report = tacopig.tests.JenkinsReport('Tacopig');



    %% Covariance function testing

    fprintf('Testing Covariance Functions______________________________________\n\n');


    covSuite = tacopig.tests.JenkinsSuite('Covariance Functions');
    report.addsuite(covSuite);
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
        thistest = tacopig.tests.JenkinsTest(current,'InitClass');


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
        thistest = tacopig.tests.JenkinsTest(current,'CheckEvalXInputDims');
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
        thistest = tacopig.tests.JenkinsTest(current,'CheckEvalVerifiesHyperLength');
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
        thistest = tacopig.tests.JenkinsTest(current,'CheckKEvalVerifiesHyperLength');
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
        thistest = tacopig.tests.JenkinsTest(current,'CheckpointvalVerifiesHyperLength');
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
        thistest = tacopig.tests.JenkinsTest(current,'CheckKevalOutputDims');
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
            thistest = tacopig.tests.JenkinsTest(current,'RegressionKevalOutput');
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
                    if (tacopig.tests.approxequal(K1,oldresult)>eqthresh)
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

        thistest = tacopig.tests.JenkinsTest(current,'KevalOutputPosDef');
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

        thistest = tacopig.tests.JenkinsTest(current,'KevalcompareEval');
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
                if (tacopig.tests.approxequal(K1,K2)>eqthresh)
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
        thistest = tacopig.tests.JenkinsTest(current,'ComparePointevalWithEval');
        if (~good_init)
            thistest.skip = 1;
        else
            try
                target = diag(Kpp)';
                obtained = covobj.pointval(X,pars);
                if (tacopig.tests.approxequal(target,obtained) > eqthresh)
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

        thistest = tacopig.tests.JenkinsTest(current,'CheckGradientsAreCorrect');
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
                        if (tacopig.tests.approxequal(target,obtained) > eqthresh)
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
    meanSuite = tacopig.tests.JenkinsSuite('Mean Functions');
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
        thistest = tacopig.tests.JenkinsTest(current,'InitClass');
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
        thistest = tacopig.tests.JenkinsTest(current,'CheckEvalVerifiesHyperLength');
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
        thistest = tacopig.tests.JenkinsTest(current,'CheckKevalOutputDims');
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
        thistest = tacopig.tests.JenkinsTest(current,'RegressionEvalOutput');
        if (~good_init)
            thistest.skip = 1;
        else
            try
                eval(sprintf('oldresult = regressions.means.%s.eval;',current));
                if (tacopig.tests.approxequal(M1,oldresult)>eqthresh)
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


        thistest = tacopig.tests.JenkinsTest(current,'CheckGradientsAreCorrect');
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
                            if (tacopig.tests.approxequal(numerical,analytic) > eqthresh)
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
    noiseSuite = tacopig.tests.JenkinsSuite('Noise Functions');
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


        thistest = tacopig.tests.JenkinsTest(current,'InitClass');
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
        thistest = tacopig.tests.JenkinsTest(current,'CheckevalOutputDims');
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
        thistest = tacopig.tests.JenkinsTest(current,'RegressionKevalOutput');
        if (~good_init)
            thistest.skip = 1;
        else
            try
                eval(sprintf('oldresult = regressions.noise.%s.eval;',current));

                if (tacopig.tests.approxequal(N1,oldresult)>eqthresh)
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
        thistest = tacopig.tests.JenkinsTest(current,'CheckGradientsAreCorrect');
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
                        if (tacopig.tests.approxequal(numerical,analytic) > eqthresh)
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


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Script to test regressor on the 1d example problem
    %  Instead of visualising, it runs a test suite and generates a report
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fprintf('Testing GP Regression framework___________________________________\n');
    % Make sure we have the optimization toolbox available
    addpath(genpath(['../../optimization']))
    %report = tacopig.tests.JenkinsReport('Tacopig');
    suite = tacopig.tests.JenkinsSuite('Test_Regressor');
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
    thistest = tacopig.tests.JenkinsTest('Regressor','Setup');
    try
        GP = tacopig.gp.Regressor;
        GP.X = X;
        GP.y = y;
        GP.MeanFn  = tacopig.meanfn.FixedMean(mean(y));
        GP.CovFn   = tacopig.covfn.SqExp();
        GP.NoiseFn = tacopig.noisefn.Stationary();
        GP.objective_function = @tacopig.objectivefn.NLML;
        GP.covpar   = 0.5*ones(1,GP.CovFn.npar(size(X,1)));
        GP.meanpar  = zeros(1,GP.MeanFn.npar(size(X,1)));
        GP.noisepar = 1e-1*ones(1,GP.NoiseFn.npar);
        GP.check();
        GP.opts.Display = 'off';
    catch e
        str = e.message;
        thistest.fail = 1;
        thistest.message = 'Failed initial setup.';
        thistest.description = str;
        fprintf('%s\n',str);
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

        thistest = tacopig.tests.JenkinsTest('Regressor','Predictive mean pre-learning.');
        if (tacopig.tests.approxequal(mf,mftarget)>eqthresh)
            thistest.fail = 1;
            thistest.message = 'Regression Test Fail';
            thistest.description = 'GP output mean (pre learning) has changed from historical value.';
            fprintf('Problem: %s\n', thistest.description);
        end
        suite.addtest(thistest);

        thistest = tacopig.tests.JenkinsTest('Regressor','Predictive variance pre-learning.');
        if (tacopig.tests.approxequal(vf,vftarget)>eqthresh)
            thistest.fail = 1;
            thistest.message = 'Regression Test Fail';
            thistest.description = 'GP output variance (pre learning) has changed from historical value.';
            fprintf('Problem: %s\n', thistest.description);
        end

        suite.addtest(thistest);

    else
        thistest = tacopig.tests.JenkinsTest('Regressor','Predictive mean pre-learning.');
        thistest.skip = 1;
        suite.addtest(thistest);
        thistest = tacopig.tests.JenkinsTest('Regressor','Predictive variance pre-learning.');
        thistest.skip = 1;
        suite.addtest(thistest);
    end


    % Pre-Learning Query
    if (good_setup)
        GP.factorisation = 'SVD';
        GP.solve(); 
        [mf, vf] = GP.query(xstar);

        thistest = tacopig.tests.JenkinsTest('Regressor','Predictive mean pre-learning with SVD factorisation.');
        if (tacopig.tests.approxequal(mf,mftarget)>eqthresh)
            thistest.fail = 1;
            thistest.message = 'Regression Test Fail';
            thistest.description = 'GP output mean (pre learning) has changed from historical value.';
            fprintf('Problem: %s\n', thistest.description);
        end
        suite.addtest(thistest);

        thistest = tacopig.tests.JenkinsTest('Regressor','Predictive variance pre-learning with SVD factorisation.');
        if (tacopig.tests.approxequal(vf,vftarget)>eqthresh)
            thistest.fail = 1;
            thistest.message = 'Regression Test Fail';
            thistest.description = 'GP output variance (pre learning) has changed from historical value.';
            fprintf('Problem: %s\n', thistest.description);
        end

        suite.addtest(thistest);

    else
        thistest = tacopig.tests.JenkinsTest('Regressor','Predictive mean pre-learning.');
        thistest.skip = 1;
        suite.addtest(thistest);
        thistest = tacopig.tests.JenkinsTest('Regressor','Predictive variance pre-learning.');
        thistest.skip = 1;
        suite.addtest(thistest);
    end



    thistest = tacopig.tests.JenkinsTest('Regressor','Learn with NLML and minfunc.');
    if (good_setup)
        GP.factorisation = 'CHOL';
        GP.learn();
        GP.solve();
        err = tacopig.tests.approxequal([GP.covpar, GP.meanpar,GP.noisepar],[1.1119    0.7811  -0.0097]);
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

    if (good_setup)
        GP.covpar = [1.1119  0.7811];
        GP.meanpar = [];
        GP.noisepar = -0.0097;
        GP.solve();
    end
    %% Generate samples from prior and posterior

    thistest = tacopig.tests.JenkinsTest('Regressor','Sample prior with fixed seed.');
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

         if (tacopig.tests.approxequal(fstar,target)>eqthresh)
            thistest.fail = 1;
            thistest.message = 'Prior sample has unexpected realisation.';
            thistest.description = 'Prior sample has changed from historical value.';
            fprintf('Problem: %s\n', thistest.description);       
        end 
    else
        thistest.skip = 1;
    end
    suite.addtest(thistest);

    thistest = tacopig.tests.JenkinsTest('Regressor','Sample posterior with fixed seed.');
    if (good_setup)
        fstar = GP.sampleposterior(xstar, 50);
        target = ...
            [-3.3448   -3.3287   -3.1678   -2.9875   -2.7723   -2.6610   -2.5796   -2.4957   -2.3547   -2.1563   -2.0280   -2.0285   -2.1607   -2.3902   -2.5087   -2.4279   -2.1966   -1.9036...
             -1.6212   -1.4226   -1.3724   -1.3864   -1.4377   -1.5112   -1.4786   -1.3090   -0.9856   -0.5933   -0.3206   -0.2412   -0.4201   -0.7638   -1.2526   -1.7597   -2.2419   -2.5897...
             -2.8252   -2.8755   -2.8625   -2.7466   -2.5878   -2.4558   -2.3459   -2.2536   -2.1564   -2.0910   -2.0797   -2.1123   -2.1007   -2.0137   -1.7543];
        if (tacopig.tests.approxequal(fstar,target)>eqthresh)
            thistest.fail = 1;
            thistest.message = 'Posterior sample has unexpected realisation.';
            thistest.description = 'Posterior sample has changed from historical value.';
            fprintf('Problem: %s\n', thistest.description);       
        end 
    else
        thistest.skip = 1;
    end
    suite.addtest(thistest);

    thistest = tacopig.tests.JenkinsTest('Regressor','Check LML Function.');
    if (good_setup)
        theta = [GP.meanpar, GP.covpar, GP.noisepar];
        [LML,LMLg] = GP.objectfun(theta);

        if (tacopig.tests.approxequal(LML,-2.6014)>eqthresh)
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

    thistest = tacopig.tests.JenkinsTest('Regressor','Check LML gradient.');
    if (good_setup)
        if (tacopig.tests.approxequal(LMLg,[-0.0271; 0.0050; -2.3410])>eqthresh)
            thistest.fail = 1;
            thistest.message = 'Objective function gradient problem.';
            thistest.description = 'Objective function has changed from historical value.';
            fprintf('Problem: %s\n', thistest.description); 
        end
    else
       thistest.skip = 1;
    end 
    fprintf('__________________________________________________________________\n');


    % Save the regression tests results as they stand
    try
        save regression_test_cache.mat regressions; % regressions
    catch e
        fprintf('Could not save regression test values... try saving manually!');
    end
    fprintf('End of test.\n');

    % wirte the report to disk
    if (nargin>0)
        try
            fid = fopen(reportpath,'w');
            fprintf(fid, report.print);
            fclose(fid);
        catch e
            fprintf('Error writing XML report: %s\n', e.message);
        end
    end
end


% Inspects a matlab error and strips out the line number link metadata for
% use in the tacopig error report
function msg = geterr(e)
     str = e.message;
     msg = [];
     active = false; % cuts off the start and the HREF
     spacecount = 0;
     newline = sprintf('\n');
     for i=1:length(str)
         if (str(i) == '<')
             active = false;
         elseif (str(i) == '>')
             active = true;
         elseif (active)

             if (str(i)==' ')
                 spacecount = spacecount+1;
                 if (spacecount ==2)
                    msg = msg(6:end);
                    msg = sprintf('%-11s ',msg);
                 else
                     msg(end+1) = ' ';
                 end
             else
                 if (str(i) ~= newline)
                    msg(end+1) = str(i); 
                 else
                    msg = [msg, ' :: '];
                 end
             end
         end
     end
     msg(end+1)  = newline;
end

% Tacopig testing structure to check that a method throws a particular error
% type under controlled conditions.
function caught = throws(in, task, desired_exception, varargin)    
    caught = 0;
    try
       eval(task);
    catch e
       caught = strcmp(e.identifier,desired_exception);

       if (~caught)
           fprintf('Warning: we got a different error - %s', e.message);
       end

    end
end