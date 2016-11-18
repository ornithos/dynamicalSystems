function [m, P] = sarkka_rts(Y, kf_m, kf_P, A, Q)

    ms = kf_m(:,end);
    Ps = squeeze(kf_P(:,:,end));
    rts_m = zeros(size(ms,1),size(Y,2));
    rts_P = zeros(size(Ps,1),size(Ps,2),size(Y,2));
    rts_m(:,end) = ms;
    rts_P(:,:,end) = Ps;
    for k=size(kf_m,2)-1:-1:1
        mp = A*kf_m(:,k);
        Pp = A*kf_P(:,:,k)*A'+Q;
        Ck = kf_P(:,:,k)*A'/Pp; 
        ms = kf_m(:,k) + Ck*(ms - mp);
        Ps = kf_P(:,:,k) + Ck*(Ps - Pp)*Ck';
        rts_m(:,k) = ms;
        rts_P(:,:,k) = Ps;
    end
    
    m = rts_m;
    P = rts_P;
end