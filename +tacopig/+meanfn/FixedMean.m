% FixedMean
% GP mean function with no parameters. Always returns its member 'value'.
% Its value is fixed during the learning phase.

classdef FixedMean < tacopig.meanfn.MeanFunc

    % Protected member variables:
    properties
        value;  % The fixed value of the Fixed mean
    end

    properties (Constant)
        ExampleUsage = 'tacopig.meanfn.FixedMean(0)'; %Instance of class created for testing
    end
    
    
    methods(Static)
        function n = npar(~)
            % Returns the number of parameters required by FixedMean()
            % This will always be zero as it does not require and parameters
            n = 0; 
        end
    end
    
    methods
        % Constructor: default to zero mean
        function this = FixedMean(val)
            % Product covariance function constructor
            %
            % GP.MeanFn = tacopig.meanfn.FixedMean(val)
            %
            % Inputs: val = The value of the fixed mean
            
           this.value = 0;
           if nargin > 0
                if (length(val)>1)
                    error('tacopig:inputInvalidLength','FixedMean cannot be a vector of length greater than 1.')
                end
               this.value = val;
           end
        end    

        % X is usually passed in, but we ignore it with ~
        function mu = eval(this, X,GP) 
            %Evaluate the mean function at locations X 
            %
            % mu = GP.MeanFn.eval(X,GP)
            %
            % GP.MeanFn is an instantiation of the FixedMean mean function class
            % Inputs:  X = D x N Input locations
            %          GP = The GP class instance can be passed to give the mean function access to its properties
            % Outputs: mu = vector of the mean at input location (1 x N)
        
            par = this.getMeanPar(GP);
            if (numel(par)~=0)
                error('tacopig:inputInvalidLength','FixedMean has no hyperparameters!')
            end
            
            mu = this.value*ones(1,size(X,2));

        end
        
        function g = gradient(this,X,par)
            %Evaluate the gradient of the mean function at locations X with respect to the parameters
            %
            % g = GP.MeanFn.gradient(X,par)
            %
            % GP.MeanFn is an instantiation of the FixedMean mean function class
            % Inputs:  X = D x N Input locations 
            %          par = The parameters of the mean function
            % Outputs: g = the gradient of the mean function at locations X with respect to the parameters (A cell of dimensionality 1 x Number of parameters. Each element is an array of dimensionality N x N)
            %
            % For this class g = []; because it has no parameters...
            g = []; % we dont have any parameters...
        end
        
        
    end
end    


