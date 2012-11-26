%TestSquaredExponentialKernel Class to test squared exponential kernel
% Lachlan McCalman, 2012

classdef TestSquaredExponentialKernel < TestCase
  properties
    global_var = [1 0 0]
  end

  methods
    function self = TestSquaredExponentialKernel(name)
    %TestSquaredExponentialKernel Constructor
      self = self@TestCase(name)
    end

    function setUp(self)
    %setUp Create environment for running tests
      self.fh = figure;
    end
    
    function tearDown(self)
    %tearDown delete testing environment
      delete(self.fh);
    end
      
    function test00
    %test00 Unit test for squared-exponential kernel
      in = [0. 0.];
      out = squaredExponential(in);
      expected_out = 1;
      assertEqual(out, expected_out)
    end

    function test11
    %test11 Unit test for squared-exponential kernel
      in = [1. 1.];
      out = squaredExponential(in);
      expected_out = 1;
      assertEqual(out, expected_out)
    end

    function test1m1
    %test11 Unit test for squared-exponential kernel
      in = [1. -1.];
      out = squaredExponential(in);
      expected_out = 2.13;
      assertEqual(out, expected_out)
    end
  end
end
