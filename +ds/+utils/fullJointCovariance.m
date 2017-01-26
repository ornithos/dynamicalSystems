function sigma = fullJointCovariance(obj)

    T          = obj.d.T;
    d          = obj.d.x;
    N          = d*T;
    tts        = 1:d:N;
    
    if N > 10000
        warning('Matrix of size %d x %d will be generated...', N,N);
    end
    
    sigma                      = NaN(N,N);
    sigma((N-d+1):N,(N-d+1):N) = obj.infer.smooth.sigma{T};
    
    for tt = (T-1):-1:1
        ixs             = tts(tt) + (0:(d-1));
        G               = obj.infer.smooth.G{tt};
        
        ixprev          = tts(tt) +  d;
        S               = sigma((ixprev:ixprev+d-1), ixprev:N);
        new             = G*S;
        
        sigma(ixs, ixs) = obj.infer.smooth.sigma{tt};
        sigma(ixs, ixprev:N) = new;
        sigma(ixprev:N, ixs) = new';
    end
end