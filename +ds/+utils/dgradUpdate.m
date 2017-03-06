function D = dgradUpdate(obj, D, P_minus, S, K, v, m, P, k)
    if nargin < 9 || isempty(k)
        k = 1;
    end
    
    assert(k == 1, 'dgrad only implemented for k = 1 at present');
    % get shortcuts to parameters
    n         = obj.d.x;
    q         = obj.d.y;
    
    A         = obj.par.A;
    H         = obj.par.H;
    Q         = obj.par.Q;
    R         = obj.par.R;
    
    invS      = inv(S);
    
    I_m       = sparse(1:q,1:q,1,q,q);
    I_n       = sparse(1:n,1:n,1,n,n);
    
    
    % Update 1
    D.v       = - full(H * D.m_minus);  % give up on sparsity here?
    D.S       = kron(H, H) * D.P_minus;
    
    %%
    % add to phi
    D.phi     = D.phi ...
                + 0.5*(kron(v',v')*kron(invS, invS) - invS(:)')*D.S ...
                - v'*invS*D.v;
    % for kk = 1:k-1: % step ahead fc
    D.G_A     = D.phi;
    
    %%
    % Update 2
    D.K       = full(kron(invS*H, I_n)*D.P_minus) - kron(I_m, P_minus*H')*kron(invS, invS)*D.S;
    D.m       = D.m_minus  + kron(v', I_n)*D.K + K*D.v;  % <----
    tmp       = kron(K*S, I_n)*D.K;
    D.P       = D.P_minus - tmp - utils.math.vecPermTranspose(n,n)*tmp - kron(K, K)*D.S;
    
    % Update 3
    D.m_minus = kron(m', I_n) + A*D.m;
    tmp       = kron(A * P, I_n);
    D.P_minus = tmp + utils.math.vecPermTranspose(n,n)*tmp + kron(A, A)*D.P;
end
    
    