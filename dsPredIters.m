niter = 100;
d     = 4;
n     = 3;
T     = 100;

testpredahead = 10;

itersLlh  = NaN(niter, 1000);
itersPred = NaN(niter, 1000);

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
    transition(randperm(d^2,round(0.5*d^2 -1e-5))) = 0;
    fprintf('\n max eigenvalue of A is %.4f.\n', max(abs(eig(transition))));
    
    emission            = 2*rand(n,d) - 1;
    emission(randperm(d*n,round(0.5*d*n))) = 0;
    
    Q                   = utils.rand.generatePSDmatrix(d,[],1e-1);
    R                   = utils.rand.generatePSDmatrix(n,[],1e-1);
    
    % generate and score
    x0                  = [];
    dsTemp  = ds.dynamicalSystem(d, n, 'x0', 1e-5, 'evolution', transition, Q, 'emission', emission, R, 'data', T);
    dsTemp.filter([],true);
    dsTemp.smooth;
    dsTemp.save('original-params');
    origLlh(kk)            = dsTemp.infer.llh;
    origPred(kk)           = sqrt(sum(sum((dsTemp.getPredictedValues(testpredahead) - ...
                                        dsTemp.y(:,(testpredahead+1):T)).^2))./dsTemp.d.T);
    
    % initialise SSID 
    opts = struct('maxiter', 0, 'ssid', true, 'verbose', false);
    dsTemp.parameterLearningEM(opts);
    
    % EM iterations (x20 per call)
    opts = struct('maxiter', 20, 'epsilon', 1e-4, 'sampleStability', 5, ...
           'multistep', 4, 'ssid', false, 'verbose', false, 'stableVerbose', false, ...
           'annealingSchedule', Inf, 'annealingIter', 100);
    
    ppb = utils.base.objProgressBar;
    ppb.newProgressBar('EM iterations: ' , 30, 'showElapsed', true);
    for ii = 1:500

        dsTemp.parameterLearningEM(opts);
        dsTemp.filter([],true);
        dsTemp.smooth;
    
        itersLlh(kk,ii)        = dsTemp.infer.llh;
        itersPred(kk,ii)       = sqrt(sum(sum((dsTemp.getPredictedValues(testpredahead) - ...
                                        dsTemp.y(:,(testpredahead+1):T)).^2))/dsTemp.d.T);
        ppb.print(ii/500);
        if ii>1 && abs(itersLlh(kk,ii)-itersLlh(kk,ii-1)) < 1e-4
            break
        end
    end
    ppb.finish;
end

pb.finish;


figure;
plot(1-bsxfun(@rdivide,itersLlh(:,1:200), origLlh')');
figure;
plot(bsxfun(@rdivide,itersPred(:,1:200), origPred')');


improveLlh  = 1-bsxfun(@rdivide,itersLlh(~tmp,1:200), itersLlh(~tmp,1));
improvePred = bsxfun(@rdivide,itersPred(~tmp,1:200), itersPred(~tmp,1));

