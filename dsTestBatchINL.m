% Setup globals
    cols          = zeros(3,3);
    cnums         = [2 4 5];
    for kk = 1:3; cols(kk,:) = utils.plot.varyColor2(cnums(kk)); end
    litecols      = utils.plot.colortint(cols, 0.8);
    cols          = utils.plot.colortint(litecols, 1.4);
    pforward      = 150;

    
doPlot  = true;

% get patient covariates
demogLkp  = tblDemog;
pxLkp     = ionlds.utils.getPatientLookup;
demogLkp.ID = pxLkp(:,2);
demogLkp  = sortrows(demogLkp, 'ID');

% setup spline
%[~,bspl] = ionlds.explore.bsplineNonlin('knots', 6, 'doPlot', false);
knots  = [-1000.0000  -70.6237   -8.3125   -1.5334    2.7607    8.7608   70.7608  1000.0000];
bspl   = utils.spline.bspline(3, knots);
h      = @(z, pars) pars.bspl.functionEval(pars.C * z + pars.bias, pars.eta);

% setup shorthand for different nonlinear optimisation schemes
nlOpts = struct;
nlOpts.CoodAsc   = struct('chgEtas', logical([0 1 1 1 1 1 1 1 1 0]), 'display', 'none', ...
                 'optimType', 'fminunc', 'bfgsSpline', false);
nlOpts.AllConstr = struct('chgEtas', logical([0 1 1 1 1 1 1 1 1 0]), 'display', 'none', ...
                 'optimType', 'fmincon', 'bfgsSpline', true, 'small_iter', 200);
nlOpts.fixEta    = struct('optimType', 'fminunc', 'fixEta', true, 'bfgsSpline', false);
nlOpts.fixC      = struct('optimType', 'fminunc', 'fixC', true, 'bfgsSpline', false);
nlOpts.fixAll    = struct('optimType', 'fminunc', 'fixC', true, 'fixEta', true, 'bfgsSpline', false);

% longest / shortest effect types
% (half life --> base):    2^(-1/n)
% (base --> half-life):    -1/(log2(a))
hlfLifeRng       = [10, 70];
hlfLifeRng       = 2.^(-1./hlfLifeRng);


%%
A          = dsBatch.par.A;
os         = dsBatch.d.y;
ss         = 4;   %dsBatch.d.x;
Q          = dsBatch.par.Q;
R          = dsBatch.par.R;

% initialise nonlinearity
etaUB      = [140, 80, 95]';
eta        = repmat([0,0,(1 + exp(-0.5*linspace(-5,5,6))).^(-1),1,1],os,1);
eta        = bsxfun(@times, eta, etaUB);
hpar       = struct('eta', eta, 'C', dsBatch.par.H, 'bias', dsBatch.par.c, 'bspl', bspl);
opts       = struct('warnings', true);

x0mu       = dsBatch.par.x0.mu;

% ____ Get IONLDS object __________________________________________
dsINLBat   = ds.ionldsBatch(ss, os, 'x0', x0mu, eye(ss)*1e-6, 'evolution', A, Q, ...
                'emission', h, hpar, R, 'data', Ys, 'control', Us, ...
                dsBatch.par.B, false, opts);

dsINLBat.smooth('ukf', [], 'forceFilter', true);
dsINLBat.save('initial-inference');

%     dsINL.useSavedParameters('initial-inference', false);  % restart
%%


% ___ Parameter Estimation _________________________________________
nlOpts                = utils.struct.addFieldToSubStructs(nlOpts, 'monotoneWarning', false);
nlOpts                = utils.struct.addFieldToSubStructs(nlOpts, 'etaUB', etaUB);
doNL                  = nlOpts.CoodAsc;

%     doNL.dbgVisualisation = true;
emOpts  = struct('maxiter', 10, 'epsilon', 5e-4, 'multistep', 4, 'ssid', false, ...
          'verbose', true, 'stableVerbose', false, 'annealingSchedule', Inf, 'annealingIter', 100, ...
          'annealingMin', 1e-3, 'diagA', true, 'diagQ', false, 'diagAconstraints', hlfLifeRng, ...
          'strictNegativeCheck', false, 'filterType', 'ukf', 'fixBias2', false, 'nonlinOptimOpts', doNL);
