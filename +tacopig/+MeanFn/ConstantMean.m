% Matrix GP Mean Function
% This is a constant mean fixed pre-learning, not a stationary mean (which would have
% a parameter)

classdef GP_ConstantMean < GP_MeanFunc

    % Protected member variables:
    properties
        value;
    end

    methods(Static)
        function n = npar(~) 
            n = 0; 
        end
    end
    
    methods
        % Constructor: default to zero mean
        function this = GP_ConstantMean(val)
           this.value = 0;
           if nargin > 0
               this.value = val';
           end
        end    

        % X is usually passed in, but we ignore it with ~
        function mu = eval(this, X,~) 
            val = this.value;
            if (numel(val) ==1)
                mu = this.value*ones(1,size(X,2));
            else
                mu = val;
            end
        end
        
        function g = gradient(this, X, par)
            g = []; % we dont have any parameters...
        end
        
        
    end
end    


