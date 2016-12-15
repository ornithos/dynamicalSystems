out = NaN(40,16);
out(:,1) = -1;
%    1     2    3    4   5   6    7    8    9   10   11    12    13     14      15     16
% {v, b}, +1, +10, +20, +1, +10, +20, +1, +10, +20, eig1, eig2, eig3, d(A,I), d(A,U), theta

plotdims       = utils.plot.subplotdims(3);
figure('Position', [250, 250, 1200, 750])

for ii = 1:40
R      = 0.2;
ss     = 3;
os     = size(Ys{ii},1);
h      = @(z, pars) ds.utilIONLDS.gen_sigmoid(pars.C * z, pars.eta);
hpar   = struct('eta', repmat([0,100,0.05,1], os, 1), 'C', 0.01*ones(os, ss));
dsINL  = ds.ionlds(ss, os, 'x0', pars{ii}.x0, pars{ii}.V0, 'evolution', [],1, 'emission', h, hpar, R, 'data', Ys{ii}, 'control', Us{ii}, true, false);

% dsINL.ssid(30);    % initialise transition matrix (not implemented SSID for control (neither has Sidiqqi). But use this paper.
% dsINL.par.B = ones(dsINL.d.x, dsINL.d.u);
% dsINL.par.Q = eye(dsINL.d.x);

dsINL.par.A = pars{ii}.A;
dsINL.par.B = pars{ii}.B;
dsINL.par.C = pars{ii}.C;
dsINL.par.Q = pars{ii}.Q;
dsINL.par.R = pars{ii}.R;
dsINL.par.emiNLParams.eta = [pars{ii}.h_param(:,[1,2,4]), zeros(os, 1)];   % add bias
dsINL.par.emiNLParams.C   = pars{ii}.h_param(:,5:end);
dsINL.filter('ukf');
dsINL.smooth('ukf');
dsINL.save('initial-inference');

% 
% if false
% % IONLDS parameter learning
% dsINL.useSavedParameters('initial-inference');
% 
% % analytic (poor precision) gradient optimisation
% dsINL.parameterLearningEM(struct('optimType', 'analytic', 'diagQ', true, 'diagR', true));
% 
% dsINL.smooth('ukf');
% dsINL.getFittedValues;
% dsINL.save('em-analytic200');
% 
% % automatic (slow) gradient optimisation
% dsINL.parameterLearningEM(struct('optimType', 'auto', 'maxiter', 25));
% dsINL.smooth('ukf');
% dsINL.getFittedValues;
% dsINL.save('em-analytic200auto25');
% end

% IONLDS check vs KG
dsINL.useSavedParameters('initial-inference', false);
dsINL.parameterLearningEM(struct('optimType', 'auto', 'maxiter', 10, 'diagQ', true, 'diagR', true));

dsINL.smooth('ukf');
dsINL.getFittedValues;
dsINL.save('em-auto10');

predictWindows = [1, 10, 20];    %-1 = smoothed, 0 = filtered, n = n-step ahead prediction
nWindows       = numel(predictWindows);

for pp = 1:nWindows
    wwSize   = predictWindows(pp);
    y        = dsINL.y;
    subplot(plotdims(1),plotdims(2),pp)
    plot(y');
    if wwSize < 0
        dsINL.getFittedValues;  
        yhat  = dsINL.yhat;
        txt   = 'smoothed';
    else
        yhat  = dsINL.getPredictedValues(wwSize);
        yhat  = [NaN(size(yhat,1), max(0,wwSize)), yhat];
        txt   = sprintf('predicted-values (+%d)', wwSize);
    end
    smse      = nansum((y - yhat).^2, 2)./(size(yhat, 2).* var(y(:,max(0,wwSize)+1:end), [], 2));
    hold on; plot(yhat'); hold off;
    title(sprintf('%s -- SMSE: %.2f%%', txt, mean(smse)*100));
    
    if size(smse,1) == 2; smse = [smse; NaN(1, size(smse,2))]; end
    out(ii, 1 + pp + [0,3,6]) = smse;
end

out(ii,11:16) = ds.utilIONLDS.getTransitionStatistics(dsINL.par.A);
% pause
end

modelframe = [nanmean(out(:,[2,5,8]),2), nanmean(out(:,[2,5,8]+1),2), nanmean(out(:,[2,5,8]+2),2), ones(40,1), out(:,11:16)];
