function out = gen_sigmoid(X, eta)

    % eta must be (1 x 4) matrix
    %         or  (n x 4) matrix
    %  matrix of parameters, where n is the dimension of x.
    
    m     = eta(:,1);
    M     = eta(:,2);
    nu    = eta(:,3);
    gamma = eta(:,4);
    
    %out   = m + (M-m)./((1 + exp(-gamma.*X)).^(1/nu));
    denom  = 1 + exp(-bsxfun(@times, gamma, X));
    denom  = bsxfun(@power, denom, 1./nu);
    
    out    = bsxfun(@rdivide, M - m, denom);
    out    = bsxfun(@plus, m, out);
    
end