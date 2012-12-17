% A standard test on all existing covariance functions

% Detects all covariance function classes and tests the following:
%   *  parse the code with matlab
%   *  instanciate the class with either no arguments or using the code 
%      in the constant string member .teststring
%   *  Check that the npar method returns a value - uses this to make a
%      random test problem with a fixed seed reset so its repeatable

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


clear all
clear functions
clc


  % Find all the covariance functions
  slash = pwd; slash = slash(1);
  cov_functions = ['+tacopig', slash, '+covfn',slash,'*.m'];
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
      
      % Create instances of these classes
      current = cov_names{i};
      
      fprintf('Testing %s:', current);
      pass = true;
      
      
      argcount = 0;
      try
         eval(['argcount = nargin(@tacopig.covfn.', current, ');']);
      catch e
         str = geterr(e);
         fprintf('\n%s\n', str);
         fprintf('Cannot continue testing %s\n',current);
         continue
      end
      currentFull = ['tacopig.covfn.', current];
      if (argcount ~=0)
          try
            eval(['teststring = ',currentFull,'.teststring;']);
            if (length(teststring) ==0)
                error('failblog');
            end
          catch e
            fprintf(['\nNo Test String for ',current,'\n']);
            fprintf('Cannot continue testing %s\n',current);
            continue;
          end
            
          fullteststring = ['covobj = tacopig.covfn.',teststring,';'];
          try
             eval(fullteststring);
          catch e
             str = geterr(e);
             fprintf('\n%s >> %s\n',fullteststring,e.message);
             fprintf('Cannot continue testing %s\n',current);
             continue
          end
      else
         try
            eval(['covobj = ', currentFull,'();']); 
         catch e
             fprintf('\n%s\nCannot continue testing %s\n',e.message,current);
             continue;
         end
      end
      
        % Output sizes (could detect an accidental scalar output etc)
        rand('seed', 50);
        npoint = 15;
        ndims = 3;
        X  = rand(ndims,npoint);
        X2 = rand(ndims,npoint-1);
        X3 = rand(ndims-1,npoint);
        npar = 0;
        try
            npar = covobj.npar(ndims);
        catch e
            fprintf('\n%s\nFurther %s tests aborted\n\n', e.message, current);
            continue;
        end
        
        % Look for standard exceptions
        pars = rand(1,npar);
        pars2 = rand(1,npar-1);
        if ~throws(covobj, 'k = in.eval(varargin{1}, varargin{2}, varargin{3});', 'tacopig:dimMismatch', X, X3, pars)
            if (pass)
                fprintf('\n');
            end
            fprintf('%s.eval is not detecting input dimensionality mismatch between x1 and x2 and should throw tacopig:dimMismatch\n', current);
            pass = false;
        end
        
        % checking number of hyperparameters
        if ~throws(covobj, 'k = in.eval(varargin{1}, varargin{2}, varargin{3});', 'tacopig:inputInvalidLength', X, X2, [pars,1])
            if (pass)
                fprintf('\n');
            end
            fprintf('%s.eval is not detecting the wrong number of hypers and should throw tacopig:inputInvalidLength\n', current);
            pass = false;
        end     
        
        if ~throws(covobj, 'k = in.Keval(varargin{1}, varargin{2});', 'tacopig:inputInvalidLength', X, pars(1:end-1))
            if (pass)
                fprintf('\n');
            end
            fprintf('%s.Keval is not detecting the wrong number of hypers and should throw tacopig:inputInvalidLength\n', current);
            pass = false;
        end     
        
        if ~throws(covobj, 'k = in.pointval(varargin{1}, varargin{2});', 'tacopig:inputInvalidLength', X, pars(1:end-1))
            if (pass)
                fprintf('\n');
            end
            fprintf('%s.pointval is not detecting the wrong number of hypers and should throw tacopig:inputInvalidLength\n', current);
            pass = false;
        end     
        
        
        % check that Keval works and is the correct size
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
            fprintf('%s >> %s\n', current, e.message);
            pass = false;
        end
        
        % and gives the same result as before
        try
            eval(sprintf('oldresult = regressions.cov.%s.Keval;',current));
            
            if (approxequal(K1,oldresult)>eqthresh)
                if (pass)
                     fprintf('\n');
                end
                fprintf('K1 doesnt match the cached regression data!\n');
                
                reply = input('Are you sure it is correct NOW? Y/N [Y]: ', 's');
                if isempty(reply)
                    reply = 'N';
                end
                if strcmpi(reply,'y')
                    eval(sprintf('regressions.cov.%s.eval = M1;',current));
                else
                    pass = false;
                end
                
                
                pass = false;
                return
                %eval(sprintf('regressions.cov.%s.Keval = K1',current));
            end
        catch e
            fprintf('No previous result for %s.Keval - saving for next time...\n', current);
            eval(sprintf('regressions.cov.%s.Keval = K1;',current));
        end
        
        
        
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
        end

        try
		    K2 = covobj.eval(X,X,pars);
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
        end				
				
        try
		    Kpp = covobj.eval(X,X,pars);
            s = size(K2);
            if ~((s(1)==npoint)&&(s(2)==npoint))
                error('%s.Keval returned a matrix of the wrong size!',current);
            end
            
            if (approxequal(K1,K2) > eqthresh)
                error(sprintf('%s.Keval(X,par) and %s.eval(X,X,par) are not giving the same answer!\n',current, current));
            end
            
        catch e
            if (pass)
                fprintf('\n');
                pass = false;
            end
            fprintf('%s\n',e.message);
		end	
        
        % check pointeval = diag
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
            fprintf('%s\n',e.message);
        end
        
        % If no gradients then warn
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
        
        % If has gradients then check that they actually correspond to
        % working gradients
        if (testgrads)
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
            end   
        end
        
        
        % passed everything?
        if (pass)
            fprintf('Passed!\n');
        end
        
        % If not gradients - issue a warning here
        
        fprintf('\n');
  end
  
    % Save the regression tests results as they stand
    try
        save regression_test_cache.mat regressions; % regressions
    catch e
        fprintf('Could not save regression test values... try saving manually!');
    end
    
  % opentoline is handy
  % containers.Map is also handy
