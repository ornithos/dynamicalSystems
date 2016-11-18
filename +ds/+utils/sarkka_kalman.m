function [m, P] = sarkka_kalman(Y, A, Q, H, R, m0, P0)

    m = m0;
    P = P0;
    kf_m = zeros(size(m,1),size(Y,2));
    kf_P = zeros(size(P,1),size(P,2),size(Y,2));
    for k=1:size(Y,2)
        m = A*m;
        P = A*P*A' + Q;
        
        S = H*P*H' + R;
        K = P*H'/S;
        m = m + K*(Y(:,k) - H*m);
        P = P - K*S*K';
        
        kf_m(:,k) = m;
        kf_P(:,:,k) = P;
    end
    
    m = kf_m;
    P = kf_P;
end