% Compares two outputs returning a positive score that increases with 
% discrepancy. A heuristic based on correlation and mean difference. 
% Usage:
%       score = approxequal(K1,K2)
function score = approxequal(K1,K2)
    if ((size(K1,1) ~=size(K2,1))||(size(K1,2) ~=size(K2,2)))
        score = inf; % not comparable
        return
    end

    mu = mean(K1(:));
    a = K1(:) - mu;
    b = K2(:) - mu;
    corr = mean(a.*b)./mean(a.*a);
    score = max( abs(1-corr), abs(mean(K1(:)) - mean(K2(:))) );
