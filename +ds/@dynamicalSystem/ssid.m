function ssid(obj, L)
    
    % dimension of input vectors
    d       = obj.d.y;
    n       = obj.d.x;
    rows    = obj.d.T - L + 1;
    Y       = zeros(d .* rows, L);
    
    % Create Hankel matrix
    for ll = 1:L
        idx     = ((ll-1)*d + 1):((ll-1+rows)*d);
        Y(:,ll) = obj.y(idx);
    end
    
    % recover the gamma matrix
    [u, s, v] = svd(Y, 'econ');
    
    gamma     = u(:,1:n) * sqrt(s(1:n,1:n));
    
    H         = gamma(1:d,:);
    gamma1    = gamma(1:end-d,:);
    gamma2    = gamma(d+1:end,:);
    A         = gamma1 \ gamma2;
    

    % Save results in object
    obj.par.H     = H;
    obj.par.A     = A;
    
end