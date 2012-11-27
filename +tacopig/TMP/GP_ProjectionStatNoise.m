% Stationary iid noise for projected sensors where the 
% size of Y is not related to the size of X.

classdef GP_ProjectionStatNoise < GP_NoiseFunc
   
    properties
        counts; 
    end
    
    % Most noise functions will be static
    methods
        function this = GP_ProjectionStatNoise(counts)
            this.counts = counts;
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
