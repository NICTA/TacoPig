classdef Stationary < tacopig.noisefn.NoiseFunc
   
    % Most noise functions will be static
    % The value of the noise is constant and learnt during training.

    methods(Static) 
        
        function n_theta = npar(~)
            % Returns the number of parameters required by the class
            % n = GP.NoiseFn.npar()
            % GP.NoiseFn is an instantiation of the Stationary noise function class
            % Always returns 1.
            n_theta = 1; % normal iid noise
        end
        
        function noise = eval(X, GP)
            % Returns the value of the noise at the location X
            %
            % noise = GP.NoiseFn.eval(X, GP)
            %
            % Gp.NoiseFn is an instantiation of the Stationary noise function class
            % Inputs: X = the input locations (D x N) 
            %         GP = The GP class instance can be passed to give the noise function access to its properties
            % Outputs : noise matrix (N x N)
            par = tacopig.noisefn.NoiseFunc.getNoisePar(GP);
            [D,N] = size(X); %N number of points in X
            noise = par^2*eye(N);
        end
        
        function [g] = gradient(X, GP)
            %Evaluate the gradient of the noise matrix at locations X with respect to the parameters
            %
            % g = GP.NoiseFn.gradient(X,GP)
            %
            % Gp.NoiseFn is an instantiation of the Stationary noise function class
            % Inputs:  X = D x N Input locations 
            % Outputs: g = the gradient of the noise function at locations X with respect to the parameters (A cell of dimensionality 1 x Number of parameters. Each element is an array of dimensionality N x N)
            %
            % For this class g is a 1 x NumberOfNoiseParameters cell array with the element being a N x N matrix (the gradient of the noise matrix with respect to the ith parameter).
         
            par = tacopig.noisefn.NoiseFunc.getNoisePar(GP);
            [D,N] = size(X); %N number of points in X
            g = {2*par*eye(N)};
        end

    end
end
