function out = dbgFreeEnergy(obj)
    T      = obj.d.T;
    d      = obj.d.x;
    n      = obj.d.y;
    
    A      = obj.par.A;
    B      = obj.par.B;
    H      = obj.par.H;
    Q      = obj.par.Q;
    R      = obj.par.R;
    
    y      = obj.y;
    u      = obj.u;
    m      = [obj.infer.smooth.x0.mu, obj.infer.smooth.mu];
    PP     = obj.infer.smooth.sigma;
    P      = zeros(d,d,T+1);
    P(:,:,1)   = obj.infer.smooth.x0.sigma;
    for ii = 2:T+1
        P(:,:,ii) = PP{ii-1};
    end
    
    GG     = obj.infer.smooth.G;
    G      = zeros(d,d,T);
    PGT    = zeros(d,d,T); 
    G(:,:,1)   = obj.infer.smooth.x0.G;
    PGT(:,:,1) = P(:,:,2) * G(:,:,1)';
    for ii = 2:T
        G(:,:,ii)   = GG{ii-1};
        PGT(:,:,ii) = P(:,:,ii+1) * G(:,:,ii)';
    end
    
    out    = zeros(5,1);
    out(1) = -0.5*T*log(det(2*pi*Q));
    out(2) = -0.5*T*log(det(2*pi*R));
    
    tmp3 =   m(:,2:T+1)*m(:,2:T+1)' + sum(P(:,:,2:T+1),3) ...
             - (m(:,2:T+1)*m(:,1:T)' + sum(PGT, 3))*A' ...
             - ((m(:,2:T+1)*m(:,1:T)' + sum(PGT, 3))*A')' ...
             + A *( m(:,1:T)*m(:,1:T)' + sum(P(:,:,1:T),3) )* A' ...
             - (m(:,2:T+1) * u') * B' ...
             - ((m(:,2:T+1) * u') * B')' ...
             + A * (m(:,1:T) * u') * B' ...
             + (A * (m(:,1:T) * u') * B')' ...
             + B * (u * u') * B';
    fprintf('Original tmp3 = \n'); disp(tmp3);
    tmp3   = 0;
    for tt = 2:T+1
        tmp3 = tmp3 + m(:,tt)*m(:,tt)' + P(:,:,tt) - (m(:,tt)*m(:,tt-1)' + P(:,:,tt)*G(:,:,tt-1)')*A' ...
               - A*(m(:,tt-1)*m(:,tt)' + G(:,:,tt-1)*P(:,:,tt)) + A*(m(:,tt-1)*m(:,tt-1)' + P(:,:,tt-1)) * A' ...
               - m(:,tt)*u(:,tt-1)'*B' - B*u(:,tt-1)*m(:,tt)' + A*(m(:,tt-1)*u(:,tt-1)')*B' ...
               + B*(u(:,tt-1)*m(:,tt-1)')*A' + B*u(:,tt-1)*u(:,tt-1)'*B';
    end
    fprintf('New tmp3 = \n'); disp(tmp3);
    
    out(3) = -0.5*trace(inv(Q) * tmp3);
    tmp4   = y*y' - H * (m(:,2:T+1)*y') - (y*m(:,2:T+1)')*H' + H*(m(:,2:T+1)*m(:,2:T+1)' + sum(P(:,:,2:T+1),3))*H';
    out(4) = -0.5*trace(inv(R) * tmp4);
    
    fullcov = ds.utils.fullJointCovariance(obj);
    detArg  = 2*pi*fullcov;
    L       = chol(detArg);
    out(5)  = 0.5*2*sum(log(diag(L))) + obj.d.T*obj.d.x./2;  % logdet = 2*sum(log(diag(chol(.))))
end