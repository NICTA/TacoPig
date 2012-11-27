%TestSquaredExponentialKernel Class to test squared exponential kernel
% Lachlan McCalman, 2012

classdef TestSqExpCovFn < TestCase
  properties
    %global_var = [1 0 0]
  end

  methods
    function self = TestSqExpCovFn(name)
    %TestSquaredExponentialKernel Constructor
      self = self@TestCase(name);
    end

    function setUp(self)
    %setUp Create environment for running tests
    end
    
    function tearDown(self)
    %tearDown delete testing environment
    end

    function testNpar1(self)
    %testNpar1 npar(1) = 2
      assertEqual(tacopig.CovFn.SqExp.npar(1), 2);
    end
    function testNpar3(self)
    %testNpar3  npar(3) = 4
      assertEqual(tacopig.CovFn.SqExp.npar(3), 4);
    end
    function testNpar0(self)
    %testNpar0  npar(0) is out of range
      f = @() tacopig.CovFn.SqExp.npar(0);
      assertExceptionThrown(f, 'tacopig:inputOutOfRange');
    end
    function testNparM1(self)
    %testNparM1 npar(-1) is out of range
      f = @() tacopig.CovFn.SqExp.npar(-1);
      assertExceptionThrown(f, 'tacopig:inputOutOfRange');
    end
    function testNpar1p5(self)
    %testNpar1p5 npar(1.5) float is invalid input 
      f = @() tacopig.CovFn.SqExp.npar(1.5);
      assertExceptionThrown(f, 'tacopig:inputInvalidType');
    end
    
    function testEval11(self)
    %testEval11 should return par(2)^2
      X1 =  1;
      X2 =  1;
      par = [1 2];
      result = tacopig.CovFn.SqExp.eval(X1, X2, par);
      answer =  4;
      assertEqual(result, answer);
    end

  end
end
