% Multi-task squared exponential for a GP model
% Assuming length scales are along axis dimensions and not affine.
% With different hyperparameters for each task. This function is 
% Called pairwise between tasks...
% Based on Melkumyan & Ramos, "Multi Kernel Gaussian Processes"
% Alistair Reid Sept 2012

classdef GP_MTSqExpCov < GP_MTCovFunc
    
    % Keep in mind, this covariance is PAIRWISE...
    % Task signal variances are now handled by the multi-task GP
    % and need not be encoded here...
    
    methods(Static) 
        
        function n_theta = npar(D)
            n_theta = D; % Per task per dimension
                         % Signal variance done externally
        end
                
        function K = eval(X1, X2, par1, par2)
            [D,N1] = size(X1);  %   number of points in X1
            N2 = size(X2,2);    %   number of points in X2
            if D~=size(X2,1)
                error('Dimensionality of X1 and X2 must be the same');
            end
            if (length(par1)~=D)
                error('Wrong number of hyperparameters 1 for MTSqExp (it doesnt use sig var)!');
            end
            if (length(par2)~=D)
                error('Wrong number of hyperparameters 2 for MTSqExp (it doesnt use sig var)!');
            end
            
            %Compute weighted squared distances:
            w = (par1.^2 + par2.^2).^-1; % effective length scales
            
            % Calculate gain between terms
            % normalisation cross length scale in the diagonal case prod is equiv to det
            gain = 2^(D/2) * prod(sqrt(abs(par1.*par2)))/sqrt(prod(par1.^2 + par2.^2));
            
            XX1 = sum(w(:,ones(1,N1)).*X1.*X1,1);
            XX2 = sum(w(:,ones(1,N2)).*X2.*X2,1);
            X1X2 = (w(:,ones(1,N1)).*X1)'*X2;
            XX1T = XX1';
            % numerical effects can drive z slightly negative 
            z = max(0,XX1T(:,ones(1,N2)) + XX2(ones(1,N1),:) - 2*X1X2);
            K = gain*exp(-0.5*z); 
        end
        
        
        function [g1, g2] = gradient(X1, X2, par1, par2)
            [D,N1] = size(X1);  %   number of points in X1
            N2 = size(X2,2);    %   number of points in X2
            if D~=size(X2,1)
                error('Dimensionality of X1 and X2 must be the same');
            end
            if (length(par1)~=D)
                error('Wrong number of hyperparameters{1} for MTSqExp (it doesnt use sig var)!');
            end
            if (length(par2)~=D)
                error('Wrong number of hyperparameters{2} for MTSqExp (it doesnt use sig var)!');
            end
            
            % D hypers per cell - no noise, no sigvar
            g1 = cell(1,D);
            g2 = cell(1,D);
            
            % Just return numerical gradients for now... to test the framework
            eval = @GP_MTSqExpCov.eval;
            eps = 1e-6;
            for i=1:D
                par1a = par1; par1a(i) = par1a(i) - eps;
                par1b = par1; par1b(i) = par1b(i) + eps;
                g1{i} = (eval(X1, X2, par1b, par2) - eval(X1, X2, par1a, par2))/(2*eps);
                
                par2a = par2; par2a(i) = par2a(i) - eps;
                par2b = par2; par2b(i) = par2b(i) + eps;
                g2{i} = (eval(X1, X2, par1, par2b) - eval(X1, X2, par1, par2a))/(2*eps);
            end
            
            return
            
            % If you need analytic gradients for this bit, can you please
            % implement: here are some notes
                % gain = 2^(D/2) * prod(sqrt(abs(par1.*par2)))/sqrt(prod(par1.^2 + par2.^2));
                % 
                % w = (par1.^2 + par2.^2).^-1; 
                % 
                % % effective length scale is...
                % LS eff sqrt(par1.^2 + par2.^2)
                %  
                %  
                % % squared exponential gradient wrt LS of each dimension
                % Kg = GP_SqExpCov.eval(X, X, par);
                % 
                % [d,n] = size(X);
                % g = cell(1,d+1);
                % for i=1:d
                %     %Compute weighted squared distance
                %     row = X(i,:);
                %     XX = row.*row;
                %     XTX = row'*row;
                %     XX = XX(ones(1,n),:);
                %     z = max(0,XX+XX'-2*XTX);
                %     w = par(i)^(-3);
                %     g{i} = w*Kg.*z;
                % end
                % g{d+1} = Kg*(2/par(end));
            

        end
       
        % Also overload the point covariance kx*x* - its trivial
        function v = pointval(x_star, par)
            %D = size(x_star,1);
            %gain = 2^(D/2) * prod(sqrt(abs(par1.*par2)))/sqrt(prod(par1.^2 + par2.^2));
            % of course gain is one because pointval implies the same
            % function!
            
            v = ones(1,size(x_star,2));% *gain
        end
        
    end
end
