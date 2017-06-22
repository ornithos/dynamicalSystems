function ssid(obj, L)
    assert(sum(sum(isnan(obj.y)))==0, 'SSID will fail with missing values at present');
    
    if ~any(obj.hasControl)
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


        % for now...
        obj.par.Q     = eye(n);
        obj.par.R     = eye(d);
    
    else
        % Using Peter Van Overschee's code; I'm too lazy to code this up
        % myself at present. It's computationally expensive to do
        % the control matrices without considerable thought.
        [A,B,H,C,K,R] = vanOverscheeSSID(obj.y,obj.u,L,obj.d.x);

        obj.par.A     = A;
        obj.par.H     = H;
        if obj.hasControl(1)
            obj.par.B = B;
        end
        if obj.hasControl(2)
            obj.par.C = C;
        end
        obj.par.R = R;
        obj.par.Q = K*R*K';
        
%         % Create control matrix if reqd.
%         if any(obj.hasControl)
%             du = obj.d.u;
%             U  = zeros(du * rows, L);
%             for ll = 1:L
%                 idx = ((ll-1)*du + 1):((ll-1+rows)*du);
%                 U(:,ll) = obj.u(:,idx);
%             end
%             projOrthU = eye(L) - U'*pinv(U*U')*U;
%         else
%             projOrthU = eye(L);
%         end
%     % Control matrices
%     if any(hasControl)
%         % No, I can't be bothered...
%         vecY = obj.y(:);
%     end

    end
end