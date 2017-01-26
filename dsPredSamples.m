niter = 60;
d     = 4;
n     = 3;
T     = 100;

testpredahead = 10;

results = zeros(niter, 6);
pb = utils.base.objProgressBar;
pb.newProgressBar('sample: ' , 30, 'showElapsed', true);

for kk = 1:niter 
    fprintf('\n');
    pb.print(kk/niter);
    fprintf('\n');
    pb.currOutputLen = 0;
    
    transition          = 2*rand(d) - 1;
    [u,s,v]             = svd(transition);
    s                   = diag(min(diag(s), 1));
    transition          = u*s*v';
    transition(randperm(d^2,round(0.5*d^2))) = 0;
    fprintf('\n max eigenvalue of A is %.4f.\n', max(abs(eig(transition))));
    
    emission            = 2*rand(n,d) - 1;
    emission(randperm(d*n,round(0.5*d*n))) = 0;
    
    Q                   = utils.rand.generatePSDmatrix(d,[],1e-1);
    R                   = utils.rand.generatePSDmatrix(n,[],1e-1);
    
    x0                  = [];
    dsTemp  = ds.dynamicalSystem(d, n, 'x0', 1e-5, 'evolution', transition, Q, 'emission', emission, R, 'data', T);
    dsTemp.filter([],true);
    dsTemp.smooth;
    dsTemp.save('original-params');

    llh.orig            = dsTemp.infer.llh;
    pred.orig           = sqrt(sum(sum((dsTemp.getPredictedValues(testpredahead) - ...
                                        dsTemp.y(:,(testpredahead+1):T)).^2)));
    
    opts = struct('maxiter', 100, 'epsilon', 1e-3, 'sampleStability', 10, ...
           'multistep', 4, 'ssid', true, 'verbose', true, 'stableVerbose', false, ...
           'annealingSchedule', Inf, 'annealingIter', 100);

    dsTemp.parameterLearningEM(opts);
    dsTemp.filter([],true);
    dsTemp.smooth;
    dsTemp.save('learned-all-EM1');
    
    llh.underfit        = dsTemp.infer.llh;
    pred.underfit       = sqrt(sum(sum((dsTemp.getPredictedValues(testpredahead) - ...
                                        dsTemp.y(:,(testpredahead+1):T)).^2)));
                                    
    opts = struct('maxiter', 1000, 'epsilon', 1e-3, 'sampleStability', 10, ...
           'multistep', 4, 'ssid', true, 'verbose', true, 'stableVerbose', false, ...
           'annealingSchedule', 0.5, 'annealingIter', 80);
    
    dsTemp.parameterLearningEM(opts);
    dsTemp.filter([],true);
    dsTemp.smooth;
    dsTemp.save('learned-all-EM2');
    
    
    llh.converge        = dsTemp.infer.llh;
    pred.converge       = sqrt(sum(sum((dsTemp.getPredictedValues(testpredahead) - ...
                                        dsTemp.y(:,(testpredahead+1):T)).^2)));
    results(kk,1:6) = [llh.orig, llh.underfit, llh.converge, ...
                       pred.orig, pred.underfit, pred.converge];
                                    
    fprintf('\n');
end

pb.finish;