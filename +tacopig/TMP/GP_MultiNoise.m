% iid noise where the input counts is used to assign different noise level
% to different points.

classdef GP_MultiNoise < GP_NoiseFunc
   
    properties
        counts; 
    end
    
    % Most noise functions will be static
    methods
        function this = GP_MultiNoise(counts)
            this.counts = counts; % Vector where each element specifies the number of observations subject to each noise parameter
        end
        
        function n_theta = npar(this,~)
            n_theta = length(this.counts);
        end
        
        function noise = eval(this, GP, par)
            N = size(GP.y,2);
            counts = this.counts;
            if (sum(counts)~=N)
                error('The noise function does not agree with the size of y');
            end
            total = sum(counts);
            ntask = length(counts);
            upto = 1;
            
            D = zeros(1,total);
            for i=1:ntask
                D(upto:upto+counts(i)-1) = par(i).^2;
                upto = upto + counts(i);
            end
            noise = diag(D);
        end
        
        function [g] = gradient(this, GP, par)
            N = size(GP.y,2);
            counts = this.counts;
            if (sum(counts)~=N)
                error('The noise function does not agree with the size of y');
            end
            total = sum(counts);
            ntask = length(counts);
            upto = 1;
            g = cell(1,ntask);
            for i=1:ntask
                D = zeros(1,total);
                D(upto:upto+counts(i)-1) = 2*par(i);
                upto = upto + counts(i);
                g{i} = diag(D);
            end
        end
    end
end
