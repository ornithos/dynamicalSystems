function obj = ssid(obj, L)
    
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
    
    % thinking about a more 'robust' method for calculating H. Turns out
    % this is way less stable due to the high power of A in these calcs.
%     [uu, ss, vv] = svd(A, 'econ');
%     ssinv     = 1./diag(ss);
%     stripGam  = zeros((rows-1)*d, n);
%     avgH = H;
%     for jj = 1:rows-1
%         stripGam((jj-1)*d+1:jj*d, :) = gamma2((jj-1)*d+1:jj*d,:) * vv * diag(ssinv.^jj) * uu';
%         avgH = H + stripGam((jj-1)*d+1:jj*d, :);
%     end
%     H = H./rows;

    % Save results in object
    obj.par.H     = H;
    obj.par.A     = A;
    
end