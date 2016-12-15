function out = gen_sigmoid(X, eta)

    % eta must be (1 x 4) matrix
    %         or  (n x 4) matrix
    %  matrix of parameters, where n is the dimension of x.
    
    m     = eta(:,1);
    M     = eta(:,2);
    nu    = eta(:,3);
    gamma = eta(:,4);
    b     = eta(:,5);
    
    %out   = m + (M-m)./((1 + exp(-gamma.*X)).^(1/nu));
    gX     = bsxfun(@times, gamma, X);
    denom  = 1 + exp(-bsxfun(@plus, gX, b));
    denom  = bsxfun(@power, denom, 1./nu);
    
    out    = bsxfun(@rdivide, M - m, denom);
    out    = bsxfun(@plus, m, out);
    
end