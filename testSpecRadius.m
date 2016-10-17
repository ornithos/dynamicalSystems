n = 100000;
d = 4;
rho = zeros(n,1);
out = zeros(n,12);
for ii = 1:n
    A = rand(d,d);
    %A = transition;
    P = utils.rand.generatePSDmatrix(d,[],true);
    X = P * A';
    
    Q = eye(d); %utils.rand.generatePSDmatrix(d,[],0.5);
    
    G = X / (A*X + Q);
    [~,uu,~] = svd(G);
    rho(ii) = max(abs(diag(uu)));
    
    out(ii,1:4) = eig(P);
    out(ii,5:8) = eig(A);
    out(ii,9)   = sum(sum(P.^2));
    out(ii,10)  = cond(P);
    out(ii,11)  = sum(sum(A.^2));
    out(ii,12)  = cond(A);
end


%%
% Find stationary covariance matrix
P = zeros(d,d);
while true
    tmp = P * obj.H';
    
    P = dsNewton.A
