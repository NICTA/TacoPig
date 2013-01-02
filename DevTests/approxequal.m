function score = approxequal(K1,K2)

    if ((size(K1,1) ~=size(K2,1))||(size(K1,2) ~=size(K2,2)))
        equal = false;
        return
    end

    mu = mean(K1(:));
    a = K1(:) - mu;
    b = K2(:) - mu;
    corr = mean(a.*b)./mean(a.*a);
    score = abs(1-corr);
