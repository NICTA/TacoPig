% Matrix GP Mean Function
% This is a constant mean fixed pre-learning, not a stationary mean (which would have
% a parameter)

classdef GP_ProjectedMean < GP_MeanFunc

    % member variables:
    properties
        sensorvalue; % what you subtract from the sensors
        cellvalue;   % what you subtract from the cells
    end

    methods(Static)
        function n = npar(~) 
            n = 0;
        end
    end
    
    methods

        function this = GP_ProjectedMean(G, Cells)
           this.sensorvalue = (G*Cells)';
           this.cellvalue = Cells;
        end    

        % X is usually passed in, but we ignore it with ~
        function mu = eval(this, ~,~) 
            mu = this.cellvalue;
        end
        
        function mu = eval_y(this, ~, ~)
            mu = this.sensorvalue;
        end
        
        function g = gradient(this, ~, ~)
            g = []; % we dont have any parameters...
        end
        
    end
end    


