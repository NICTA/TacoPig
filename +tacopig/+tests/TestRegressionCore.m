% TestSquaredExponentialKernel Class to test squared exponential kernel
% Lachlan McCalman, 2012

classdef TestRegressionCore < TestCase

  % What do we want to test today?
  % 
  % set verbose to false
    
    
  properties
    %global_var = [1 0 0]
  end

  methods
      
    function self = TestRegressionCore(name)
    %TestSquaredExponentialKernel Constructor
      self = self@TestCase(name);
    end

    function setUp(self)
        %setUp Create environment for running tests
        fprintf('setUp\n');
    end
    
    function tearDown(self)
        %tearDown delete testing environment
        fprintf('tearDown\n');
    end

    function testOne(self)
    %testNpar1 npar(1) = 2
      fprintf('testOne\n');
      %assertEqual(tacopig.CovFn.SqExp.npar(1), 2);
    end
    
    function testTwo(self)
        fprintf('testTwo\n');
    end
    
    
    %testNpar0  npar(0) is out of range
      %f = @() tacopig.CovFn.SqExp.npar(0);
      %assertExceptionThrown(f, 'tacopig:inputOutOfRange');
  end
  
end
