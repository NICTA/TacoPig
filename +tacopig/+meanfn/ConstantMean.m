% ConstantMean
% GP mean function with no parameters. Always returns its member 'value'.
classdef ConstantMean < tacopig.meanfn.MeanFunc

    % Protected member variables:
    properties
        value;
    end

    properties (Constant)
        ExampleUsage = 'tacopig.meanfn.ConstantMean(0)';
    end
    
    
    methods(Static)
        function n = npar(~) 
            n = 0; 
        end
    end
    
    methods
        % Constructor: default to zero mean
        function this = ConstantMean(val)
           this.value = 0;
           if nargin > 0
               this.value = val';
           end
        end    

        % X is usually passed in, but we ignore it with ~
        function mu = eval(this, X,GP) 
            par = this.getMeanPar(GP);
            if (numel(par)~=0)
                error('tacopig:inputInvalidLength','ConstantMean has no hyperparameters!')
            end
            
            val = this.value;
            if (numel(val) ==1)
                mu = this.value*ones(1,size(X,2));
            else
                mu = val;
            end
        end
        
        function g = gradient(this,X,par)
            g = []; % we dont have any parameters...
        end
        
        
    end
end    


