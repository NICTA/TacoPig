classdef GP_StatNoise < GP_NoiseFunc
   
    % Most noise functions will be static
    methods(Static) 
        
        function n_theta = npar(~)
            n_theta = 1; % normal iid noise
        end
        
        function noise = eval(in1, par)
            if isa(in1, 'GP_STD')
                X = in1.X;
            else
                X = in1;
            end
            [D,N] = size(X); %N number of points in X
            noise = par^2*eye(N);
        end
        
        function [g] = gradient(in1, par)
            if isa(in1, 'GP_STD')
                X = in1.X;
            else
                X = in1;
            end
            [D,N] = size(X); %N number of points in X
            g = {2*par*eye(N)};
        end

    end
end
