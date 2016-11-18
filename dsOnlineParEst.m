%%
% Online Experiments
% progress = utils.base.objProgressBar;
% progress.newProgressBar('Simulation in progress: ', 30, true, true);
emOpts   = struct('maxiter', 2000, 'epsilon', 1e-4, 'dbg', false);
numSim   = 10;
miterNums= 200; [100,200,500];
N        = 5;
D        = 6;

% N, D, orig llh, batch llh, batch iter, batch time, online (100) - [llh,niter,time] ,
% online (500) - [llh,niter,time], online (1000) - [llh,niter,time] ;
results  = NaN(numSim, 15);
online   = cell(0,3);
% 
% rng(1200);
for ii = 1:numSim
    
    A = 2*rand(N)-1;
    [U,S,V] = svd(A, 'econ');
    S       = max(-0.8,min(0.8, S));
    A       = U*S*V';
    A(randperm(N^2, floor(N^2/2))) = 0;
    
    q       = rand(1);
    Q       = q*utils.rand.generatePSDmatrix(N) + (1-q)*eye(N);
    
    H       = rand(D, N);
    H(randperm(N*D, floor(D*N/2))) = 0;
    
    r       = rand(1);
    R       = r*utils.rand.generatePSDmatrix(D) + (1-r)*eye(D);
    
    dsTemp  = ds.dynamicalSystem(N, D, 'x0', 1e-5, 'evolution', A, Q, 'emission', H, R, 'data', 140);
    dsTemp  = dsTemp.filter('',true);
    dsTemp  = dsTemp.smooth;
    dsTemp  = dsTemp.save('original-params');
    
    dsTemp.par.A = eye(N);
    dsTemp.par.Q = eye(N);
    dsTemp.par.H = rand(D, N) - 1;
    
    if false
        tmpTic                = tic;
        [dsEM, llh, niter]    = dsTemp.parameterLearningEM(emOpts);
        results(ii,1:5)       = [N, D, dsTemp.getSaved('original-params').infer.llh, llh(end), niter];
        results(ii,6)         = toc(tmpTic);
    end
    
    % online optimisation
    for jj = 1:numel(miterNums)
        onlineOpts            = struct('maxiter', miterNums(jj), 'epsilon', 1e-3, 'dbg', false);
        [~,onlineHist]        = dsTemp.parameterLearnOnline([], onlineOpts);
        results(ii,(7:9)+(jj-1)*3)       = [onlineHist(end,1), sum(onlineHist(:,2)), onlineHist(end,3)];
        online{ii,1}          = onlineHist;
    end
%     progress.print(ii/numSim);
    fprintf('------ Completed %4d of %d simulations ---------\n', ii, numSim);
end
% progress.finish;

